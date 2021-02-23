#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <qcd.h>

#define REUNIT_INTERVAL 20

double gauge_fix_action(qcd_gaugeField *u,qcd_int_4 OrT)
{
	
	qcd_int_8 i;
	qcd_int_4 dir;
	qcd_complex_16 gf_action = (qcd_complex_16){ 0.0,0.0 };
	qcd_complex_16 gf_action_sum = (qcd_complex_16){ 0.0,0.0 };
	qcd_int_4 ndir = 4 - OrT; // 判断求和方向的个数，用于做归一化
	
	for(i=0;i<u->geo->lV;i++)
		for (dir = OrT; dir < 4; dir++)
		{
			gf_action.re += u->D[i][dir][0][0].re + u->D[i][dir][1][1].re + u->D[i][dir][2][2].re;
			gf_action.im += u->D[i][dir][0][0].im + u->D[i][dir][1][1].im + u->D[i][dir][2][2].im;
		}

	MPI_Allreduce(&gf_action.re, &gf_action_sum.re, 2, MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
//	if (u->geo->myid == 0) printf("succeed gauge_fix_action %e %e\n",gf_action_sum.re/ (3 * ndir * u->geo->V),gf_action_sum.im/ (3 * ndir * u->geo->V));

	return gf_action_sum.re/(3*ndir*u->geo->V);

}

void printmatrix3x3(qcd_complex_16 m[3][3], char* w)
{
	int i, j;
	printf(w);
	printf("\n");
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			printf("%e %e    ", m[i][j].re, m[i][j].im);
		}
		printf("\n");
	}
}





void create_rotation(qcd_gaugeTransformation * G, qcd_gaugeField* u,qcd_int_4 OrT) //G(x)=sum_mu u_mu(x) + u_mu(x-mu)^daggar
{
	
	qcd_int_8 i;
	qcd_int_4 dir;

	qcd_complex_16 tmp[3][3];
	qcd_zeroGaugeTransformation(G);
	qcd_communicateGaugePM(u);
	qcd_waitall(u->geo);

	for (i = 0; i < u->geo->lV; i++)
		for (dir = OrT; dir < 4; dir++)
		{
			qcd_addAdjoint3x3(tmp, u->D[i][dir], u->D[u->geo->minus[i][dir]][dir]);
			qcd_add3x3(G->D[i], G->D[i], tmp);
		}
}


void create_rotationsu2(qcd_gaugeTransformation* Gsu2, qcd_gaugeTransformation* G,\
	qcd_int_4 p,qcd_int_4 q,qcd_int_4 parity, qcd_real_8 overrelax)
{
	
	qcd_int_8 i,j;
	qcd_real_8 a0, a1, a2, a3, asq, a0sq;
	qcd_real_8 x,xdr,r;
	qcd_int_4 lt, lx, ly, lz;

	qcd_zeroGaugeTransformation(Gsu2);

	for (lt = 0; lt < G->geo->lL[0]; lt++)
	for (lx = 0; lx < G->geo->lL[1]; lx++)
	for (ly = 0; ly < G->geo->lL[2]; ly++)
	for (lz = 0; lz < G->geo->lL[3]; lz++)
	{
		i = qcd_LEXIC(lt, lx, ly, lz, G->geo->lL);
		if ((lt + lx + ly + lz) % 2 == parity)
		{
			a0 = G->D[i][p][p].re + G->D[i][q][q].re;
			a1 = G->D[i][p][q].im + G->D[i][q][p].im;
			a2 = G->D[i][p][q].re - G->D[i][q][p].re;
			a3 = G->D[i][p][p].im - G->D[i][q][q].im;

			asq = a1 * a1 + a2 * a2 + a3 * a3;
			a0sq = a0 * a0;
			x = overrelax;
			r = sqrt(a0sq + x * x * asq);
			xdr = x / r;
			a0 = a0 / r; a1 = a1 * xdr; a2 = a2 * xdr; a3 = a3 * xdr;

			Gsu2->D[i][p][p] = (qcd_complex_16){ a0,-a3 };
			Gsu2->D[i][p][q] = (qcd_complex_16){ -a2,-a1 };
			Gsu2->D[i][q][p] = (qcd_complex_16){ a2,-a1 };
			Gsu2->D[i][q][q] = (qcd_complex_16){ a0,a3 };

			for (j = 0; j < 3; j++)
			{
				if (j != p && j != q)
				{
					Gsu2->D[i][j][j] = (qcd_complex_16){ 1,0 };
				}
			}
		}
		else
		{
			qcd_unit3x3(Gsu2->D[i]);
		}
		
	
	}
}


void rotate(qcd_gaugeTransformation* Gsu2, qcd_gaugeField* u)// U_mu(x)=G(x) U_mu(x) G^daggar(x+mu)
{	
	qcd_int_8 i;
	qcd_int_4 dir;
	qcd_complex_16 tmp[3][3];
	 qcd_communicateTransformationM(Gsu2);
	 qcd_waitall(Gsu2->geo);
	 for(i=0;i<u->geo->lV;i++)
		 for (dir = 0; dir < 4; dir++)
		 {
			
			 qcd_MUL3x3(tmp, Gsu2->D[i], u->D[i][dir]);
			 qcd_MULADJOINT3x3(u->D[i][dir], tmp, Gsu2->D[u->geo->plus[i][dir]]);		
		 }
}

void gaugefixstep(qcd_gaugeField* u,qcd_real_8 *current_av,qcd_int_4 OrT, qcd_real_8 overrelax)
{
	
	qcd_gaugeTransformation G;
	qcd_gaugeTransformation Gsu2;
	qcd_int_8 i;
	qcd_int_4 p, q,parity;
	qcd_real_8 x, xdr, r;
	qcd_complex_16 su2[2][2];

	qcd_initGaugeTransformation(&G,u->geo);
	qcd_initGaugeTransformation(&Gsu2,u->geo);

	for (parity = 0; parity <= 1; parity++)
	{
		
		create_rotation(&G, u,OrT);

		create_rotationsu2(&Gsu2, &G, 0, 1,parity,overrelax);
		rotate(&Gsu2, u);
	
		create_rotationsu2(&Gsu2, &G, 1, 2,parity, overrelax);
		rotate(&Gsu2, u);
		
		create_rotationsu2(&Gsu2, &G, 0, 2,parity, overrelax);
		rotate(&Gsu2, u);		

	}
	*current_av = gauge_fix_action(u, OrT);
		
	qcd_destroyGaugeTransformation(&G);
	qcd_destroyGaugeTransformation(&Gsu2);
}

void reunitarize_su3(qcd_complex_16 g[3][3])
{
	qcd_real_8 ar;
	qcd_complex_16 a01;
	qcd_int_4 i;

	// 正规化第0行
	ar = qcd_NORMSQUARED(g[0][0]) + qcd_NORMSQUARED(g[0][1]) + qcd_NORMSQUARED(g[0][2]);
	ar = 1.0 / sqrt(ar);
	g[0][0] = qcd_CSCALE(g[0][0], ar);
	g[0][1] = qcd_CSCALE(g[0][1], ar);
	g[0][2] = qcd_CSCALE(g[0][2], ar);

	// 第一行正交第0 行
	a01 = qcd_CADD(qcd_CADD(qcd_CADJOINTMUL(g[0][0], g[1][0]), qcd_CADJOINTMUL(g[0][1], g[1][1]))
		, qcd_CADJOINTMUL(g[0][2], g[1][2]));   //u_0^dagger * u_1

	for (i = 0; i < 3; i++)
	{
		g[1][i] = qcd_CSUB(g[1][i], qcd_CMUL(g[0][i], a01));
	} // u_1=u_1- (u_0^dagger *u_1 ) u_0

	//正规化第一行
	ar = qcd_NORMSQUARED(g[1][0]) + qcd_NORMSQUARED(g[1][1]) + qcd_NORMSQUARED(g[1][2]);
	ar = 1.0 / sqrt(ar);
	g[1][0] = qcd_CSCALE(g[1][0], ar);
	g[1][1] = qcd_CSCALE(g[1][1], ar);
	g[1][2] = qcd_CSCALE(g[1][2], ar);

	// 第二行是第一行和第二行的叉乘再共轭 
	g[2][0] = qcd_CSUB(qcd_CMUL(g[0][1], g[1][2]), qcd_CMUL(g[0][2], g[1][1]));
	g[2][0] = qcd_CONJ(g[2][0]);
	g[2][1] = qcd_CSUB(qcd_CMUL(g[0][2], g[1][0]), qcd_CMUL(g[0][0], g[1][2]));
	g[2][1] = qcd_CONJ(g[2][1]);
	g[2][2] = qcd_CSUB(qcd_CMUL(g[0][0], g[1][1]), qcd_CMUL(g[0][1], g[1][0]));
	g[2][2] = qcd_CONJ(g[2][2]);

}//end reunitarize_su3

void  reunitarize(qcd_gaugeField* u)
{
	
	qcd_int_8 i;
	qcd_int_4 dir;
	for (i = 0; i < u->geo->lV; i++)
		for (dir = 0; dir < 4; dir++)
		{
			reunitarize_su3(u->D[i][dir]);
		}
	
} //end reunitarize


void ContinuumField(qcd_gaugeField* A, qcd_gaugeField* u, qcd_int_4 OrT)
//G(x)=sum_mu u_mu(x) - u_mu(x)^daggar then traceless
{

	qcd_int_8 i;
	qcd_int_4 dir;

	qcd_complex_16 tmp[3][3];
	qcd_complex_16 tr;
	qcd_complex_16 unit[3][3];


	for (i = 0; i < u->geo->lV; i++)
	{
		for (dir = OrT; dir < 4; dir++)
		{
			qcd_subAdjoint3x3(A->D[i][dir], u->D[i][dir], u->D[i][dir]);
			//if (u->geo->myid == 0 &&dir==0) printmatrix3x3(A->D[0][0], "A.D[0][0]");

			tr = qcd_trace3x3(A->D[i][dir]);
			//if (u->geo->myid == 0 && dir == 0) printf("trace re %e im %e \n", tr.re, tr.im);
			tr.im = tr.im / 3; //tr 实部为零
			//if (u->geo->myid == 0 && dir == 0) printf("trace re %e im %e \n", tr.re, tr.im);
			qcd_unit3x3(unit);
			qcd_cScale3x3(unit, tr);
			//if (u->geo->myid == 0 && dir == 0) printmatrix3x3(unit, "Unit");
			qcd_sub3x3(A->D[i][dir], A->D[i][dir], unit);
			//if (u->geo->myid == 0 && dir == 0) printmatrix3x3(A->D[0][0], "A.D[0][0]");
		}

	}
}


qcd_real_8 gaugecheck(qcd_gaugeField* u, qcd_int_4 OrT)
{
	qcd_gaugeField A;
	qcd_int_8 i;
	qcd_int_4 dir, j, k;
	qcd_complex_16 tmp[3][3];
	qcd_real_8 checkzero = 0;
	qcd_real_8 checkzero2;
	qcd_int_4 ndir = 4 - OrT;
	qcd_initGaugeField(&A,u->geo);

	ContinuumField(&A, u, OrT);
	qcd_communicateGaugePM(&A);
	qcd_waitall(A.geo);
	
	//if (u->geo->myid == 0) printmatrix3x3(A.D[0][0], "A.D[0][0]");
	for (i = 0; i < u->geo->lV; i++)
	{
		qcd_zero3x3(tmp);
		for (dir = OrT; dir < 4; dir++)
			for (j = 0; j < 3; j++)
				for (k = 0; k < 3; k++)
				{
					tmp[j][k] =qcd_CADD(tmp[j][k], qcd_CSUB(A.D[i][dir][j][k], A.D[A.geo->minus[i][dir]][dir][j][k]));
				}

		for (j = 0; j < 3; j++)
			for (k = 0; k < 3; k++)
			{
				checkzero +=  qcd_NORM(tmp[j][k]);
			}
		
	}
	
	MPI_Allreduce(&checkzero, &checkzero2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	qcd_destroyGaugeField(&A);
	return checkzero2/A.geo->V/ndir;
}



/*gauge type 0 or 1 : Landau or Coulomb */
void gaugefix(int gauge_type, qcd_gaugeField* u, qcd_real_8 tol, qcd_int_4 max_iter, qcd_real_8 overrelax)
{
	qcd_int_4 gauge_iter;
	qcd_real_8 current_av = 0.0, old_av = 0.0, del_av = 0.0;
	qcd_int_4 OrT = gauge_type; // 0 包括T方向Landau ，1 不包括T方向，Coulomb
	qcd_real_8 plaq;
	qcd_real_8 checkzero; 
	for (gauge_iter = 0; gauge_iter < max_iter; gauge_iter++)
	{
		gaugefixstep(u,&current_av,OrT,overrelax);

		if ((gauge_iter % REUNIT_INTERVAL) == (REUNIT_INTERVAL - 1))
		{
			reunitarize( u);
		}

		if (gauge_iter != 0)
		{
			del_av = current_av - old_av;
			if (fabs(del_av) < tol) break;
		}
		old_av = current_av;

		checkzero = gaugecheck(u, OrT);

		if (u->geo->myid == 0) printf("step %d gf action %.8e delta %.8e,the gauge check %e  \n",
			gauge_iter, current_av, del_av, checkzero);

		//plaq = qcd_calculatePlaquette(u);
		//if (u->geo->myid == 0) printf("plaq %e \n", plaq);
		
	}

	reunitarize(u);
}


