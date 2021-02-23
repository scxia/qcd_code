#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <qcd.h>
#include <qcd_projectors.h>


void PropagatorFourierTransform3D(qcd_complex_16 res2[4][4][3][3], qcd_propagator* u, qcd_int_4 mom[3], qcd_int_4 x_src[4], qcd_int_4 lt)
{
	qcd_int_4 lx, ly, lz, x, y, z;
	qcd_int_4 mu, nu, c1, c2;
	qcd_int_8 v;
	qcd_real_8 tmp;
	qcd_complex_16 C2;
	qcd_complex_16 res[4][4][3][3];

	memset(&(res[0][0][0][0].re), 0, 4 * 4 * 3 * 3 * sizeof(qcd_complex_16));

	if (lt >= 0 && lt < u->geo->lL[0])
	{
		for (lx = 0; lx < u->geo->lL[1]; lx++)
		for (ly = 0; ly < u->geo->lL[2]; ly++)
		for (lz = 0; lz < u->geo->lL[3]; lz++)
		{

			v = qcd_LEXIC(lt, lx, ly, lz, u->geo->lL);
			x = lx + u->geo->Pos[1] * u->geo->lL[1] - x_src[1];
			y = ly + u->geo->Pos[2] * u->geo->lL[2] - x_src[2];
			z = lz + u->geo->Pos[3] * u->geo->lL[3] - x_src[3];
			tmp = (((double)mom[0] * x) / u->geo->L[1] + ((double)mom[1] * y) / u->geo->L[2] + ((double)mom[2] * z) / u->geo->L[3]) * 2 * M_PI;
			C2 = (qcd_complex_16){ cos(tmp), -sin(tmp) }; //e^{-ipx}
			for (mu = 0; mu < 4; mu++)
			for (nu = 0; nu < 4; nu++)
			for (c1 = 0; c1 < 3; c1++)
			for (c2 = 0; c2 < 3; c2++)
			{
				res[mu][nu][c1][c2] = qcd_CADD(res[mu][nu][c1][c2], qcd_CMUL(u->D[v][mu][nu][c1][c2], C2));
			}

		}
	}
	MPI_Allreduce(&res[0][0][0][0].re, &res2[0][0][0][0].re, 4 * 4 * 3 * 3 * 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

}

void PropagatorFourierTransform4D(qcd_complex_16 res2[4][4][3][3], qcd_propagator* u, qcd_int_4 mom[4], qcd_int_4 x_src[4])
{
	qcd_int_4 lx, ly, lz,lt, x, y, z,t;
	qcd_int_4 mu, nu, c1, c2;
	qcd_int_8 v;
	qcd_real_8 tmp;
	qcd_complex_16 C2;
	qcd_complex_16 res[4][4][3][3];

	memset(&(res[0][0][0][0].re), 0, 4 * 4 * 3 * 3 * sizeof(qcd_complex_16));

	for (lt = 0; lt < u->geo->lL[0]; lt++)
	for (lx = 0; lx < u->geo->lL[1]; lx++)
	for (ly = 0; ly < u->geo->lL[2]; ly++)
	for (lz = 0; lz < u->geo->lL[3]; lz++)
	{

		v = qcd_LEXIC(lt, lx, ly, lz, u->geo->lL);
		t = lt + u->geo->Pos[0] * u->geo->lL[0] - x_src[0];
		x = lx + u->geo->Pos[1] * u->geo->lL[1] - x_src[1];
		y = ly + u->geo->Pos[2] * u->geo->lL[2] - x_src[2];
		z = lz + u->geo->Pos[3] * u->geo->lL[3] - x_src[3];
		tmp = (((double) mom[0] * t) /u->geo->L[0]+ ((double)mom[1] * x) / u->geo->L[1] + 
			((double)mom[2] * y) / u->geo->L[2] + ((double)mom[3] * z) / u->geo->L[3]) * 2 * M_PI;
		C2 = (qcd_complex_16){ cos(tmp), -sin(tmp) }; //e^{-ipx}
		for (mu = 0; mu < 4; mu++)
		for (nu = 0; nu < 4; nu++)
		for (c1 = 0; c1 < 3; c1++)
		for (c2 = 0; c2 < 3; c2++)
		{
			res[mu][nu][c1][c2] = qcd_CADD(res[mu][nu][c1][c2], qcd_CMUL(u->D[v][mu][nu][c1][c2], C2));
		}

	}
	
	MPI_Allreduce(&res[0][0][0][0].re, &res2[0][0][0][0].re, 4 * 4 * 3 * 3 * 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

}


void BlockMatrixFourierTransform3D(qcd_complex_16 res2[4][4], qcd_geometry *geo,qcd_complex_16(*block)[4][4], qcd_int_4 mom[3], qcd_int_4 x_src[4], qcd_int_4 lt)
{
	qcd_int_4 lx, ly, lz, x, y, z;
	qcd_int_4 mu, nu, ku, lu;
	qcd_int_8 v;
	qcd_real_8 tmp;
	qcd_complex_16 C2;
	qcd_complex_16 res[4][4];

	memset(&(res[0][0].re), 0, 4 * 4  * sizeof(qcd_complex_16));

	if (lt >= 0 && lt < geo->lL[0]) //
	{
		for (lx = 0; lx < geo->lL[1]; lx++)
		for (ly = 0; ly < geo->lL[2]; ly++)
		for (lz = 0; lz < geo->lL[3]; lz++)
		{

			v = qcd_LEXIC0(lx, ly, lz, geo->lL);
			x = lx + geo->Pos[1] * geo->lL[1] - x_src[1];
			y = ly + geo->Pos[2] * geo->lL[2] - x_src[2];
			z = lz + geo->Pos[3] * geo->lL[3] - x_src[3];
			tmp = (((double)mom[0] * x) / geo->L[1] + ((double)mom[1] * y) / geo->L[2] + ((double)mom[2] * z) / geo->L[3]) * 2 * M_PI;
			C2 = (qcd_complex_16){ cos(tmp), -sin(tmp) }; //e^{-ipx}
			for (mu = 0; mu < 4; mu++)
			for (nu = 0; nu < 4; nu++)
			{
				res[mu][nu] = qcd_CADD(res[mu][nu], qcd_CMUL(block[v][mu][nu], C2));
			}

		}
	}
	MPI_Allreduce(&res[0][0].re, &res2[0][0].re, 4 * 4 * 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

}
void BlockFourierTransform3D(qcd_complex_16 *res2, qcd_geometry* geo, qcd_complex_16 *block, qcd_int_4 mom[3], qcd_int_4 x_src[4], qcd_int_4 lt)
{
	qcd_int_4 lx, ly, lz, x, y, z;
	qcd_int_4 mu, nu, ku, lu;
	qcd_int_8 v;
	qcd_real_8 tmp;
	qcd_complex_16 C2;
	qcd_complex_16 res;

	res = (qcd_complex_16){ 0,0 };

	if (lt >= 0 && lt < geo->lL[0]) //
	{
		for (lx = 0; lx < geo->lL[1]; lx++)
		for (ly = 0; ly < geo->lL[2]; ly++)
		for (lz = 0; lz < geo->lL[3]; lz++)
		{

			v = qcd_LEXIC0( lx, ly, lz, geo->lL);
			x = lx + geo->Pos[1] * geo->lL[1] - x_src[1];
			y = ly + geo->Pos[2] * geo->lL[2] - x_src[2];
			z = lz + geo->Pos[3] * geo->lL[3] - x_src[3];
			tmp = (((double)mom[0] * x) / geo->L[1] + ((double)mom[1] * y) / geo->L[2] + ((double)mom[2] * z) / geo->L[3]) * 2 * M_PI;
			C2 = (qcd_complex_16){ cos(tmp), -sin(tmp) }; //e^{-ipx}
			res = qcd_CADD(res, qcd_CMUL(block[v], C2));
						

		}
	}
	MPI_Allreduce(&res.re, res2, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

}

void PionTwoPointWallSource(qcd_propagator* u1, qcd_propagator* u2, qcd_int_4 mom1[3], qcd_int_4 mom2[3], \
	qcd_int_4 tstart, qcd_int_4 tend, qcd_int_4 src_t, FILE* fp_corr_p)
{
	qcd_int_4 t, lt;
	qcd_int_4 mu, nu, c1, c2;
	qcd_int_4 x_src[4] = { src_t,0,0,0 };
	qcd_complex_16 res;
	qcd_complex_16 res1[4][4][3][3];
	qcd_complex_16 res2[4][4][3][3];
	for (t = tstart; t < tend; t++)
	{
		lt = ((t + src_t) % u1->geo->L[0]) - u1->geo->Pos[0] * u1->geo->lL[0];
		PropagatorFourierTransform3D(res1, u1, mom1, x_src, lt);
		PropagatorFourierTransform3D(res2, u2, mom2, x_src, lt);

		res = (qcd_complex_16){ 0,0 };
		for (mu = 0; mu < 4; mu++)
		for (nu = 0; nu < 4; nu++)
		for (c1 = 0; c1 < 3; c1++)
		for (c2 = 0; c2 < 3; c2++)
		{
			res = qcd_CADD(res, qcd_CMUL(res1[mu][nu][c1][c2], qcd_CONJ(res2[mu][nu][c1][c2])));
		}

		if (u1->geo->myid == 0)
		{
			fprintf(fp_corr_p, "%d %+e %+e\n", t, res.re, res.im);
			printf("time is %d \n", t);
		}
	}
}

void PionTwoPointMomSmearing(qcd_propagator* propp, qcd_propagator* propm, qcd_geometry* geo,qcd_int_4 tstart,qcd_int_4 tend, \
	qcd_int_4 mom[3], qcd_int_4 x_src[4], FILE* fp_corr_p)
{
	qcd_int_8 v3, v ;
	qcd_int_4 lx, ly, lz, lt;
	qcd_int_4 mu, nu;
	qcd_int_4 c1,c1p;
	qcd_complex_16 *block;
	qcd_complex_16 corr;
	qcd_int_4 t;

	block = malloc(geo->lV3  * sizeof(qcd_complex_16));

	for (t = tstart; t < tend; t++)
	{
		memset(block, 0, geo->lV3  * sizeof(qcd_complex_16));

		lt = ((t + x_src[0]) % geo->L[0]) - geo->Pos[0] * geo->lL[0]; 

		if (lt >= 0 && lt < geo->lL[0])
		{

			for (lx = 0; lx < geo->lL[1]; lx++)
			for (ly = 0; ly < geo->lL[2]; ly++)
			for (lz = 0; lz < geo->lL[3]; lz++)
			{
				v = qcd_LEXIC(lt, lx, ly, lz, geo->lL);
				v3 = qcd_LEXIC0(lx, ly, lz, geo->lL);

				for (mu = 0; mu < 4; mu++)
				for (nu = 0; nu < 4; nu++)								
				for (c1 = 0; c1 < 3; c1++)
				for (c1p = 0; c1p < 3; c1p++)
				{
						block[v3] = qcd_CADD(block[v3],
							qcd_CMUL(qcd_CONJ(propm->D[v][mu][nu][c1p][c1]), propp->D[v][mu][nu][c1p][c1]));
				}
			}
		}
	//	if (geo->myid == 0) printf("%+e %+e \n", block[0].re, block[0].im);
		BlockFourierTransform3D(&corr, propm->geo, block, mom, x_src, lt);

		if (geo->myid == 0)
		{		
			fprintf(fp_corr_p, "%d %d %d %d %+e %+e\n", t,mom[0], mom[1], mom[2],  corr.re, corr.im);
		}
	}

	free(block);
}
/*
void PionTwoPointAxial(qcd_propagator* propp, qcd_propagator* propm, qcd_geometry* geo, qcd_int_4 tstart, qcd_int_4 tend, \
	qcd_int_4 mom[3], qcd_int_4 x_src[4], FILE* fp_corr_p)
{
	qcd_int_8 v3, v;
	qcd_int_4 lx, ly, lz, lt;
	qcd_int_4 mu, nu, lu, ku;
	qcd_int_4 c1, c1p;
	qcd_complex_16* block;
	qcd_complex_16 corr;
	qcd_int_4 t;

	block = malloc(geo->lV3 * sizeof(qcd_complex_16));

	for (t = tstart; t < tend; t++)
	{
		memset(block, 0, geo->lV3 * sizeof(qcd_complex_16));

		lt = ((t + x_src[0]) % geo->L[0]) - geo->Pos[0] * geo->lL[0]; //这里需要更一般化

		if (lt >= 0 && lt < geo->lL[0])
		{

			for (lx = 0; lx < geo->lL[1]; lx++)
				for (ly = 0; ly < geo->lL[2]; ly++)
					for (lz = 0; lz < geo->lL[3]; lz++)
					{
						v = qcd_LEXIC(lt, lx, ly, lz, geo->lL);
						v3 = qcd_LEXIC0(lx, ly, lz, geo->lL);

						for (mu = 0; mu < 4; mu++)
						for (nu = 0; nu < 4; nu++)
						for (c1 = 0; c1 < 3; c1++)
						for (c1p = 0; c1p < 3; c1p++)
						{
							block[v3] = qcd_CADD(block[v3],
								qcd_CMUL(qcd_CMUL( qcd_CMUL(qcd_GAMMA[0][ku][mu],qcd_CONJ(propm->D[v][mu][nu][c1p][c1])),\
									qcd_GAMMA[0][nu][lu]),propp->D[v][ku][lu][c1p][c1]));
						}
					}
		}
		if (geo->myid == 0) printf("%+e %+e \n", block[0].re, block[0].im);
		BlockFourierTransform3D(&corr, propm->geo, block, mom, x_src, lt);

		if (geo->myid == 0)
		{
			fprintf(fp_corr_p, "%d %d %d %d %+e %+e\n", t, mom[0], mom[1], mom[2], corr.re, corr.im);
		}
	}

	free(block);
}

*/


void PionFormFactorWallSource(qcd_propagator *uprop1, qcd_propagator *uprop2, qcd_propagator *dprop1,qcd_propagator *dprop2,\
	qcd_geometry* geo, qcd_int_4 t_start, qcd_int_4 t_stop, qcd_int_4 momtf[3], qcd_int_4 tsource,FILE* fp_corr_p)
{
	qcd_int_4 t, x, y, z, lt, lx, ly, lz;
	qcd_int_4 x1, y1, z1, i;
	qcd_int_4 mu, nu, mup, nup;
	qcd_int_4 al, be, alp, bep;
	qcd_int_4 a, b, ap, bp;
	qcd_int_4 bperp;
	qcd_int_8 v3, v, vp, v3p;
	qcd_complex_16 C, C2;
	qcd_real_8 tmp, tmp2;
	qcd_complex_16  corr2;
	qcd_int_4 ctr;
	qcd_int_4 test = 0;
	qcd_int_4 x_src[4] = { 0,0,0,0 };
	qcd_complex_16  *block;

	block = malloc(geo->lV3 * sizeof(qcd_complex_16) );
	//if (geo->myid == 0) printf("here1 \n");

	for (t = t_start; t <= t_stop; t++)
	{
		lt = ((t + tsource) % geo->L[0]) - geo->Pos[0] * geo->lL[0];
		for (bperp = -geo->L[1] / 2; bperp <= geo->L[1] / 2; bperp++)
		{
			for (mu = 0; mu < 4; mu++)
			for (nu = 0; nu < 4; nu++)
			for (mup = 0; mup < 4; mup++)
			for (nup = 0; nup < 4; nup++)
			{

				for (v3 = 0; v3 < geo->lV3; v3++)   //set blocks to zero
				{
					block[v3] = (qcd_complex_16){ 0, 0 };
				}

				if (lt >= 0 && lt < geo->lL[0])
				{
					for (lx = 0; lx < geo->lL[1]; lx++)
					for (ly = 0; ly < geo->lL[2]; ly++)
					for (lz = 0; lz < geo->lL[3]; lz++)
					{
						v = qcd_LEXIC(lt, lx, ly, lz, geo->lL);
						vp = qcd_LEXIC(lt, (lx + bperp + geo->L[1]) % geo->L[1], ly, lz, geo->lL);
						v3 = qcd_LEXIC0(lx, ly, lz, geo->lL);

						for (al = 0; al < 4; al++) //for(be=0;be<4;be++)																					
						for (alp = 0; alp < 4; alp++) //for (bep = 0; bep < 4; bep++)
						{
							be = al;
							bep = alp;
							for (a = 0; a < 3; a++)
							for (b = 0; b < 3; b++)
							for (ap = 0; ap < 3; ap++)
							for (bp = 0; bp < 3; bp++)
							{
							
								block[v3] = qcd_CADD(block[v3],
									qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_GAMMA[5][bep][bep], qcd_CMUL(qcd_CONJ(dprop2->D[vp][mu][bep][b][ap]), qcd_GAMMA[5][mu][mu]))
										, qcd_CMUL(uprop2->D[vp][nu][al][b][a], qcd_GAMMA[5][al][be]))
										, qcd_CMUL(qcd_CMUL(qcd_GAMMA[5][be][be], qcd_CMUL(qcd_CONJ(uprop1->D[v][mup][be][bp][a]), qcd_GAMMA[5][mup][mup]))
											, qcd_CMUL(dprop1->D[v][nup][alp][bp][ap], qcd_GAMMA[5][alp][bep]))));

							}
						}
					}

				}


				BlockFourierTransform3D(&corr2, geo, block, momtf, x_src, lt);

				if (geo->myid == 0)
				{
					fprintf(fp_corr_p, "%d %d %d %d %d %d %d %+e %+e \n", momtf[2], t, bperp, mu, nu, mup, nup, corr2.re, corr2.im);
				}
			}
		}


		
	}



}



void StapleWilsonLine(qcd_complex_16 Staple[3][3], int lt, int lx, int ly, int lz, int l, int bprep, qcd_geometry* geo, qcd_gaugeField* u)
{
	int v, vp;
	int ll, lbprep;

	qcd_complex_16 temp[3][3];

	qcd_unit3x3(Staple);
	qcd_unit3x3(temp);

	for (ll = 0; ll < l; ll++) //Z 方向
	{
		vp = qcd_LEXIC(lt, lx, ly, (lz - ll - 1 + geo->L[3]) % geo->L[3], geo->lL);
		qcd_mulAdjoint3x3(temp, Staple, u->D[vp][3]);
		qcd_copy3x3(Staple, temp);

	}
	//  bprep 为正执行 //X 方向
	for (lbprep = 0; lbprep < bprep; lbprep++)
	{
		vp = qcd_LEXIC(lt, (lx + lbprep + geo->L[1]) % geo->L[1], ly, (lz - l + geo->L[3]) % geo->L[3], geo->lL);
		qcd_mul3x3(temp, Staple, u->D[vp][1]);
		qcd_copy3x3(Staple, temp);
	}
	//bprep 为负执行
	for (lbprep = 0; lbprep > bprep; lbprep--)
	{
		vp = qcd_LEXIC(lt, (lx + lbprep - 1 + geo->L[1]) % geo->L[1], ly, (lz - l + geo->L[3]) % geo->L[3], geo->lL);
		qcd_mulAdjoint3x3(temp, Staple, u->D[vp][1]);
		qcd_copy3x3(Staple, temp);
	}

	for (ll = 0; ll < l; ll++)
	{
		vp = qcd_LEXIC(lt, (lx + bprep + geo->L[1]) % geo->L[1], ly, (lz - l + ll + geo->L[3]) % geo->L[3], geo->lL);
		qcd_mul3x3(temp, Staple, u->D[vp][3]);
		qcd_copy3x3(Staple, temp);
	}

	qcd_dagger3x3(Staple);
	//	return wilsonline;
}



void PionQuasiWaveFunctionWallSource(qcd_propagator * prop1,qcd_propagator *prop2,qcd_gaugeField *u,qcd_geometry * geo,\
	qcd_int_4 mom[3],qcd_int_4 src_t, FILE* fp_corr_p)
{
	qcd_int_8 v3, v, vp;
	qcd_int_4 lx, ly, lz, lt;
	qcd_int_4 mu, nu, ku, lu;
	qcd_int_4 c1, c2, c1p, c2p;
	qcd_complex_16 Staple[3][3];
	qcd_complex_16(* block)[4][4];
	qcd_complex_16 corr[4][4];
	qcd_int_4 l, bprep;
	qcd_int_4 x_src[4] = { src_t,0,0,0 };

	block = malloc(geo->lV3 * 4 * 4 * sizeof(qcd_complex_16));

	for (bprep = -12; bprep <= 12; bprep++)
	for (l = 0; l <= 12; l++)
	{
		 memset(block, 0, geo->lV3 * 4 * 4 * sizeof(qcd_complex_16));

		lt = ((6 + x_src[0]) % geo->L[0]) - geo->Pos[0] * geo->lL[0]; //这里需要更一般化

		if (lt >= 0 && lt < geo->lL[0])
		{

			for (lx = 0; lx < geo->lL[1]; lx++)
			for (ly = 0; ly < geo->lL[2]; ly++)
			for (lz = 0; lz < geo->lL[3]; lz++)
			{
				v = qcd_LEXIC(lt, lx, ly, lz, geo->lL);
				vp = qcd_LEXIC(lt, (lx + bprep + geo->L[1]) % geo->L[1], ly, lz, geo->lL);
				v3 = qcd_LEXIC0(lx, ly, lz, geo->lL);
				StapleWilsonLine(Staple, lt, lx, ly, lz, l, bprep, geo, u);


				for (mu = 0; mu < 4; mu++)
				for (nu = 0; nu < 4; nu++)
				{
					for (lu = 0; lu < 4; lu++)
					for (ku = 0; ku < 4; ku++)
					for (c1 = 0; c1 < 3; c1++)
					for (c1p = 0; c1p < 3; c1p++)
					for (c2p = 0; c2p < 3; c2p++)
					{


						block[v3][mu][nu] = qcd_CADD(block[v3][mu][nu],
							qcd_CMUL(qcd_CMUL(qcd_CONJ(prop2->D[vp][ku][lu][c2p][c1p]), Staple[c2p][c1]),
								qcd_CMUL(qcd_GAMMA[5][ku][mu], prop1->D[v][nu][lu][c1][c1p])));

					}
						
				}


			}
		}	
		
		BlockMatrixFourierTransform3D(corr, u->geo,block, mom, x_src, lt);
				
		if (geo->myid == 0)
		{
			for(mu=0;mu<4;mu++)
			for(nu=0;nu<4;nu++)
			fprintf(fp_corr_p, "%d %d %d %d %d %d %d %+e %+e\n", bprep,l,mom[0],mom[1],mom[2],mu,nu,corr[mu][nu].re,corr[mu][nu].im);
		}
	}

	free(block);
}

void QuarkLocalOperatorRenormalization(qcd_propagator* prop, FILE* fp_corr_p)
{
	qcd_int_4 mu, nu, ku, lu;
	qcd_int_4 c1, c2, c1p, c2p;
	qcd_int_4 lx, ly, lz, lt;
	qcd_int_8 v;
	qcd_complex_16 block,corr;

	for (c1 = 0; c1 < 3; c1++)
	for (c2 = 0; c2 < 3; c2++)
	for (mu = 0; mu < 4; mu++)	
	for (lu = 0; lu < 4; lu++)
	for (nu = 0; nu < 4; nu++)
	for (ku = 0; ku < 4; ku++)
	{
			block = (qcd_complex_16){ 0,0 };
			for (lt = 0; lt < prop->geo->lL[0]; lt++)
			for (lx = 0; lx < prop->geo->lL[1]; lx++)
			for (ly = 0; ly < prop->geo->lL[2]; ly++)
			for (lz = 0; lz < prop->geo->lL[3]; lz++)
			//for (c1 = 0; c1 < 3; c1++)
			for (c1p = 0; c1p < 3; c1p++)
			{
				v = qcd_LEXIC(lt, lx, ly, lz, prop->geo->lL);

				block = qcd_CADD(block, qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_GAMMA[5][mu][mu], qcd_CONJ(prop->D[v][nu][mu][c1p][c1])),
					qcd_GAMMA[5][nu][nu]), prop->D[v][ku][lu][c1p][c2]));
			}
			MPI_Allreduce(&block, &corr, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

			if (prop->geo->myid == 0)
			{
				fprintf(fp_corr_p, "%d %d %d %d %d %d %+e %+e \n",c1,c2, mu, lu, nu, ku, corr.re, corr.im);
			}
	}

		
}



void QuarkQuasiWaveFunctionWallSource(qcd_propagator* prop, qcd_gaugeField* u, qcd_geometry* geo, FILE* fp_corr_p)
{
	qcd_int_8 v3, v, vp;
	qcd_int_4 lx, ly, lz, lt;
	qcd_int_4 mu, nu, ku, lu;
	qcd_int_4 c1, c2, c1p, c2p;
	qcd_complex_16 Staple[3][3];
	qcd_complex_16 block[4][4][4][4];
	qcd_complex_16 corr[4][4][4][4];
	qcd_complex_16 block2, corr2;
	qcd_int_4 l, bprep;
	qcd_complex_16 blocktest;

	qcd_propagator propga5dagger;
	qcd_initPropagator(&propga5dagger, geo);
	qcd_zeroPropagator(&propga5dagger);

	
	for (v = 0; v<geo->lV;v++)
	for (mu = 0; mu < 4; mu++)
	for (ku = 0; ku < 4; ku++)
	for(c2p=0;c2p<3;c2p++)
	for(c1=0;c1<3;c1++)
	{
		propga5dagger.D[v][ku][mu][c2p][c1] = qcd_CMUL(qcd_GAMMA[5][ku][ku], qcd_CMUL(qcd_CONJ(prop->D[v][mu][ku][c2p][c1]),
			qcd_GAMMA[5][mu][mu]));
	}


	for (bprep = -12; bprep <= 12; bprep++)
	for (l = 0; l <= 0; l++)
	{	
		if (geo->myid == 0) printf("b %d L %d \n", bprep, l);
		for (c1 = 0; c1 < 3; c1++)
		{
			c2 = c1;
		
				memset(block, 0, 4 * 4 *4*4* sizeof(qcd_complex_16));
				memset(corr, 0, 4 * 4 * 4 * 4 * sizeof(qcd_complex_16));
			for (mu = 0; mu < 4; mu++)
			for (nu = 0; nu < 4; nu++)
			for (lu = 0; lu < 4; lu++)
			for (ku = 0; ku < 4; ku++)
			{
				for (lt = 0; lt < geo->lL[0]; lt++)
				for (lx = 0; lx < geo->lL[1]; lx++)
				for (ly = 0; ly < geo->lL[2]; ly++)
				for (lz = 0; lz < geo->lL[3]; lz++)
				{
					v = qcd_LEXIC(lt, lx, ly, lz, geo->lL);
					vp = qcd_LEXIC(lt, (lx + bprep + geo->L[1]) % geo->L[1], ly, lz, geo->lL);
					StapleWilsonLine(Staple, lt, lx, ly, lz, l, bprep, geo, u);
						
						for (c1p = 0; c1p < 3; c1p++)
						for (c2p = 0; c2p < 3; c2p++)
						{

							block[ku][lu][mu][nu] = qcd_CADD(block[ku][lu][mu][nu],
								qcd_CMUL(qcd_CMUL(propga5dagger.D[vp][ku][mu][c2p][c1], Staple[c2p][c1p]),prop->D[v][nu][lu][c1p][c2]));
						}

				}

			}
	
			MPI_Allreduce(&(block[0][0][0][0].re), &(corr[0][0][0][0].re), 4 * 4 * 4 * 4 * 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

				if (geo->myid == 0)
				{
					for (ku = 0; ku < 4; ku++)
					for (lu = 0; lu < 4; lu++)
					for (mu = 0; mu < 4; mu++)
					for (nu = 0; nu < 4; nu++)
						fprintf(fp_corr_p, "%d %d %d %d %d %d %d %d %+e %+e\n", bprep, l,c1,c2, ku,lu,mu, nu, corr[ku][lu][mu][nu].re, corr[ku][lu][mu][nu].im);
				}
		}
	}

	qcd_destroyPropagator(&propga5dagger);
}

void QuarkWaveFunctionFourierTransformation(qcd_propagator* prop, qcd_int_4 mom[4], FILE * fp_corr_p)
{
	qcd_complex_16 res[4][4][3][3];
	qcd_int_4 x_src[4] = { 0,0,0,0 };
	qcd_int_4 mu, nu,c1,c2;
	
	PropagatorFourierTransform4D(res, prop, mom, x_src);
	if (prop->geo->myid == 0)
	{
		for (mu = 0; mu < 4; mu++)
		for (nu = 0; nu < 4; nu++)
		for (c1 = 0; c1 < 3; c1++)
		for (c2 = 0; c2 < 3; c2++)
		{		
			fprintf(fp_corr_p, "%d %d %d %d %+e %+e \n", mu, nu,c1,c2, res[mu][nu][c1][c2].re, res[mu][nu][c1][c2].im);
		}
	}
}


