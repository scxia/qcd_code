#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <qcd.h>
#include <qcd_projectors.h>


void twop_pion_soft_function(qcd_propagator *dprop_smear, qcd_propagator *upropsink_smear, qcd_geometry *geo,qcd_int_4 x_src[4], qcd_int_4 x_sink[4], qcd_int_4 mom[3], FILE *fp_corr_p)
{
	qcd_int_4 lx1[4];
	qcd_int_4 lx2[4];
	qcd_int_4 i,mu,nu,a,b;
	qcd_complex_16 u[4][4][3][3],d[4][4][3][3];
	qcd_complex_16 u2[4][4][3][3], d2[4][4][3][3];
	qcd_complex_16 corr;
	qcd_int_8 v;
	

	memset(u, 0, 4 * 4 * 3 * 3 * sizeof(qcd_complex_16));
	memset(d, 0, 4 * 4 * 3 * 3 * sizeof(qcd_complex_16));
	memset(u2, 0, 4 * 4 * 3 * 3 * sizeof(qcd_complex_16));
	memset(d2, 0, 4 * 4 * 3 * 3 * sizeof(qcd_complex_16));

	for (i = 0;  i < 4;i++)
	{
		lx1[i] = x_src[i] - geo->lL[i] * geo->Pos[i];
		lx2[i] = x_sink[i] - geo->lL[i] * geo->Pos[i];
	}

	if (lx2[0] >= 0 && lx2[0] < geo->lL[0] && lx2[1] >= 0 && lx2[1] < geo->lL[1] && \
		lx2[2] >= 0 && lx2[2] < geo->lL[2] && lx2[3] >= 0 && lx2[3] < geo->lL[3])
	{
		v = qcd_LEXIC(lx2[0], lx2[1], lx2[2], lx2[3], geo->lL);
		for (mu = 0; mu < 4; mu++)
			for (nu = 0; nu < 4; nu++)
				for (a = 0; a < 3; a++)
					for (b = 0; b < 3; b++)
					{
						u[mu][nu][a][b] = dprop_smear->D[v][mu][nu][a][b];

					}
	}


	if (lx1[0] >= 0 && lx1[0] < geo->lL[0] && lx1[1] >= 0 && lx1[1] < geo->lL[1] && \
		lx1[2] >= 0 && lx1[2] < geo->lL[2] && lx1[3] >= 0 && lx1[3] < geo->lL[3])
	{
		v = qcd_LEXIC(lx1[0], lx1[1], lx1[2], lx1[3], geo->lL);
		for (mu = 0; mu < 4; mu++)
			for (nu = 0; nu < 4; nu++)
				for (a = 0; a < 3; a++)
					for (b = 0; b < 3; b++)
					{
						d[mu][nu][a][b] = upropsink_smear->D[v][mu][nu][a][b];

					}
	}

	MPI_Allreduce(&(u[0][0][0][0].re), &(u2[0][0][0][0].re), 4 * 4 * 3 * 3 * 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&(d[0][0][0][0].re), &(d2[0][0][0][0].re), 4 * 4 * 3 * 3 * 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);


	if (geo->myid == 0)
	{
		for (mu = 0; mu < 4; mu++)
			for (nu = 0; nu < 4; nu++)
				for (a = 0; a < 3; a++)
					for (b = 0; b < 3; b++)
					{
						corr = qcd_CMUL(qcd_CMUL(u2[mu][nu][a][b], qcd_GAMMA[5][nu][nu]), qcd_CMUL(d2[nu][mu][b][a], qcd_GAMMA[5][mu][mu]));
					}
		//这里需要添加傅里叶变换系数，现在零动量先忽略
		fprintf(fp_corr_p, "%d %e %e \n",(x_sink[0]-x_src[0]+geo->L[0])%geo->L[0] , corr.re,corr.im);

	}


}






void form_factor_soft_function(qcd_propagator *dprop, qcd_propagator *upropsink, qcd_geometry *geo, \
	qcd_int_4 t_start, qcd_int_4 t_stop, qcd_int_4 x_src[4], qcd_int_4 x_sink[4], qcd_int_4 mom[3],qcd_int_4 momtf[3], FILE *fp_corr_p)
{
	qcd_int_4 t, x, y, z, lt, lx, ly, lz;
	qcd_int_4 x1, y1, z1,i;
	qcd_int_4 mu, nu, mup, nup;
	qcd_int_4 al, be, alp, bep;
	qcd_int_4 a, b, ap, bp;
	qcd_int_4 bperp;
	qcd_int_8 v3,v,vp,v3p;
	qcd_complex_16 C,C2;
	qcd_real_8 tmp,tmp2;
	qcd_complex_16 corr[4][4][4][4],corr2[4][4][4][4];
	qcd_int_4 ctr;
	qcd_int_4 test = 0;

	qcd_complex_16 (*block)[4][4][4][4];

	block = malloc(geo->lV3 * sizeof(qcd_complex_16) * 16 * 16);
	//if (geo->myid == 0) printf("here1 \n");

	for (bperp=-geo->L[1]/2;bperp<=geo->L[1]/2;bperp++)
	{

		for (t = t_start; t <= t_stop; t++)
		{
			lt = ((t + x_src[0]) % geo->L[0]) - geo->Pos[0] * geo->lL[0];

			for (v3 = 0; v3 < geo->lV3; v3++)   //set blocks to zero
			{
				for (mu = 0; mu < 4; mu++)
					for (nu = 0; nu < 4; nu++)
						for (mup = 0; mup < 4; mup++)
							for (nup = 0; nup < 4; nup++)
							{
								block[v3][mu][nu][mup][nup] = (qcd_complex_16) { 0, 0 };
							}
			}
			//test = 0;
			//if (geo->myid == 0) printf("here2 and lt = %d \n",lt);
			if (lt >= 0 && lt < geo->lL[0])
			{
				for (lx = 0; lx < geo->lL[1]; lx++)
					for (ly = 0; ly < geo->lL[2]; ly++)
						for (lz = 0; lz < geo->lL[3]; lz++)
						{
							v= qcd_LEXIC(lt, lx, ly, lz, geo->lL);
							vp = qcd_LEXIC(lt, (lx + bperp + geo->L[1]) % geo->L[1], ly, lz, geo->lL);
							v3 = qcd_LEXIC0(lx, ly, lz, geo->lL);
							for (mu = 0; mu < 4; mu++)
								for (nu = 0; nu < 4; nu++)
									for (mup = 0; mup < 4; mup++)
										for (nup = 0; nup < 4; nup++)
										{
											for(al=0;al<4;al++) //for(be=0;be<4;be++)																					
													for(alp=0;alp<4;alp++) //for (bep = 0; bep < 4; bep++)
													{
														be = al;
														bep = alp;
															for (a = 0; a < 3; a++)
																for (b = 0; b < 3; b++)
																	for (ap = 0; ap < 3; ap++)
																		for (bp = 0; bp < 3; bp++)
																		{

																			block[v3][mu][nu][mup][nup] = qcd_CADD(block[v3][mu][nu][mup][nup],
																				qcd_CMUL(qcd_CMUL(qcd_CMUL(qcd_GAMMA[5][bep][bep], qcd_CMUL(qcd_CONJ(dprop->D[vp][mu][bep][b][ap]), qcd_GAMMA[5][mu][mu]))
																				, qcd_CMUL(upropsink->D[vp][nu][al][b][a], qcd_GAMMA[5][al][be]))
																				, qcd_CMUL(qcd_CMUL(qcd_GAMMA[5][be][be], qcd_CMUL(qcd_CONJ(upropsink->D[v][mup][be][bp][a]), qcd_GAMMA[5][mup][mup]))
																					, qcd_CMUL(dprop->D[v][nup][alp][bp][ap], qcd_GAMMA[5][alp][bep]))));
																			//test++;

																		}
													}
										}

						}



			}
			//if (geo->myid == 0) printf("here3 and test is %d\n",test);
			for (i = 0; i < 4; i++)
			{
				momtf[2] = i * 2;
				mom[2] = i;
				/*x1 = x_src[1] - x_sink[1];
				y1 = x_src[2] - x_sink[2];
				z1 = x_src[3] - x_sink[3];
				tmp2 = (((double)mom[0] * x1) / geo->L[1] + ((double)mom[1] * y1) / geo->L[2] + ((double)mom[2] * z1) / geo->L[3]) * 2 * M_PI;
				C = (qcd_complex_16) { cos(tmp2), -sin(tmp2) };
				if (geo->myid == 0) printf("%+e %+e \n", C.re, C.im); */

				for (mu = 0; mu < 4; mu++)
					for (nu = 0; nu < 4; nu++)
						for (mup = 0; mup < 4; mup++)
							for (nup = 0; nup < 4; nup++)
							{
								corr[mu][nu][mup][nup] = (qcd_complex_16) { 0, 0 };
							}
				//if (geo->myid == 0) printf("here4 \n");
				for (lx = 0; lx < geo->lL[1]; lx++)
					for (ly = 0; ly < geo->lL[2]; ly++)
						for (lz = 0; lz < geo->lL[3]; lz++)
						{
							//C = (qcd_complex_16) { 1, 0 };


							v3 = qcd_LEXIC0(lx, ly, lz, geo->lL);
							x = lx + geo->Pos[1] * geo->lL[1] - x_src[1];
							y = ly + geo->Pos[2] * geo->lL[2] -x_src[2];
							z = lz + geo->Pos[3] * geo->lL[3] -x_src[3];
							tmp = (((double)momtf[0] * x) / geo->L[1] + ((double)momtf[1] * y) / geo->L[2] + ((double)momtf[2] * z) / geo->L[3]) * 2 * M_PI;
							C2 = (qcd_complex_16) { cos(tmp), sin(tmp) };

							for (mu = 0; mu < 4; mu++)
								for (nu = 0; nu < 4; nu++)
									for (mup = 0; mup < 4; mup++)
										for (nup = 0; nup < 4; nup++)
										{
											corr[mu][nu][mup][nup] = qcd_CADD(corr[mu][nu][mup][nup], qcd_CMUL(block[v3][mu][nu][mup][nup], C2));
										}
						}

				//if (geo->myid == 0) printf("here5 \n");
				MPI_Allreduce(&(corr[0][0][0][0].re), &(corr2[0][0][0][0].re), 2 * 16 * 16, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);



				//if (geo->myid == 0) printf("here6 \n");

				if (geo->myid == 0)
				{
					x = x_sink[1] - x_src[1];
					y = x_sink[2] - x_src[2];
					z = x_sink[3] - x_src[3];
					tmp = (((double)mom[0] * x) / geo->L[1] + ((double)mom[1] * y) / geo->L[2] + ((double)mom[2] * z) / geo->L[3]) * 2 * M_PI;
					C2 = (qcd_complex_16) { cos(tmp), -sin(tmp) };
					if (geo->myid == 0) printf("%d %d %d %+e %+e \n", x,y,z,C2.re, C2.im);

					for (mu = 0; mu < 4; mu++)
						for (nu = 0; nu < 4; nu++)
							for (mup = 0; mup < 4; mup++)
								for (nup = 0; nup < 4; nup++)
								{
									corr2[mu][nu][mup][nup] = qcd_CMUL(corr2[mu][nu][mup][nup], C2);
								}
					//if (geo->myid == 0) printf("here7 \n");
					for (mu = 0; mu < 4; mu++)
						for (nu = 0; nu < 4; nu++)
						{
							//fprintf(fp_corr_p, "%d %d ", bperp, t);
							for (mup = 0; mup < 4; mup++)
								for (nup = 0; nup < 4; nup++)
								{
									fprintf(fp_corr_p, "%d %d %d %d %d %d %d %+e %+e \n", mom[2],bperp, t, mu, nu, mup, nup, corr2[mu][nu][mup][nup].re, corr2[mu][nu][mup][nup].im);
								}
							//fprintf(fp_corr_p, "\n");
						}

					//	if (geo->myid == 0) printf("here8 \n");


				}
			}


		}
	}



}

void create_wilson_line(qcd_complex_16 wilsonline[3][3],int lt,int lx, int ly,int lz, int l, int bprep,qcd_geometry *geo,qcd_gaugeField *u)
{
	int v, vp;
	int ll, lbprep;
	
	qcd_complex_16 temp[3][3];

	qcd_unit3x3(wilsonline);
	qcd_unit3x3(temp);

	for (ll = 0; ll < l; ll++) //Z 方向
	{
		vp = qcd_LEXIC(lt, lx, ly, (lz-ll -1 + geo->L[3])%geo->L[3], geo->lL);
		qcd_mulAdjoint3x3(temp, wilsonline, u->D[vp][3]);
		qcd_copy3x3(wilsonline, temp);

	}
	//  bprep 为正执行
	for (lbprep = 0; lbprep < bprep; lbprep++)
	{
		vp = qcd_LEXIC(lt, (lx + lbprep + geo->L[1]) % geo->L[1], ly, (lz-l + geo->L[3]) % geo->L[3], geo->lL);
		qcd_mul3x3(temp, wilsonline, u->D[vp][1]);
		qcd_copy3x3(wilsonline, temp);
	}
	//bprep 为负执行
	for (lbprep = 0; lbprep > bprep; lbprep--)
	{
		vp = qcd_LEXIC(lt, (lx + lbprep-1 + geo->L[1]) % geo->L[1], ly, (lz - l + geo->L[3]) % geo->L[3], geo->lL);
		qcd_mulAdjoint3x3(temp, wilsonline, u->D[vp][1]);
		qcd_copy3x3(wilsonline, temp);
	}

	for (ll = 0; ll < l; ll++)
	{
		vp= qcd_LEXIC(lt, (lx + bprep + geo->L[1]) % geo->L[1], ly, (lz-l+ ll + geo->L[3]) % geo->L[3], geo->lL);
		qcd_mul3x3(temp, wilsonline, u->D[vp][3]);
		qcd_copy3x3(wilsonline, temp);
	}

	qcd_dagger3x3(wilsonline);
//	return wilsonline;
}



qcd_complex_16 create_wilson_loop( int lt, int lx, int ly, int lz, int l, int bprep, qcd_geometry *geo, qcd_gaugeField *u)
{
	int v, vp;
	int ll, lbprep;
	qcd_complex_16 wilsonline[3][3];
	qcd_complex_16 temp[3][3];
	qcd_complex_16 wilsonloopValue,su3test;

	qcd_unit3x3(wilsonline);


	for (ll = 0; ll < l; ll++) //Z 方向
	{
		vp = qcd_LEXIC(lt, lx, ly, (lz - ll - 1 + geo->L[3]) % geo->L[3], geo->lL);
		qcd_mulAdjoint3x3(temp, wilsonline, u->D[vp][3]);
		qcd_copy3x3(wilsonline, temp);

	}
	// 暂时 bprep 只能为正
	for (lbprep = 0; lbprep < bprep; lbprep++)
	{
		vp = qcd_LEXIC(lt, (lx + lbprep + geo->L[1]) % geo->L[1], ly, (lz - l + geo->L[3]) % geo->L[3], geo->lL);
		qcd_mul3x3(temp, wilsonline, u->D[vp][1]);
		qcd_copy3x3(wilsonline, temp);
	}

	for (ll = 0; ll < l; ll++)
	{
		vp = qcd_LEXIC(lt, (lx + bprep + geo->L[1]) % geo->L[1], ly, (lz - l + ll + geo->L[3]) % geo->L[3], geo->lL);
		qcd_mul3x3(temp, wilsonline, u->D[vp][3]);
		qcd_copy3x3(wilsonline, temp);
	}


	for (lbprep = 0; lbprep < bprep; lbprep++)
	{
		vp = qcd_LEXIC(lt, (lx + bprep-1-lbprep + geo->L[1]) % geo->L[1], ly, lz , geo->lL);
		qcd_mulAdjoint3x3(temp, wilsonline, u->D[vp][1]);
		qcd_copy3x3(wilsonline, temp);
	}

	qcd_dagger3x3(wilsonline);
	wilsonloopValue = qcd_trace3x3(wilsonline);

/*	qcd_mulAdjoint3x3(temp, wilsonline, wilsonline);
	su3test=qcd_trace3x3(temp);

	printf("%e %e\n", su3test.re, su3test.im);*/

	return wilsonloopValue;
}

/*qcd_complex_16 create_wilson_loop(int lt, int lx, int ly, int lz, int l, int bprep, qcd_geometry *geo, qcd_gaugeField *u)
{
	int v, vp;
	int ll, lbprep;
	qcd_complex_16 wilsonline[3][3];
	qcd_complex_16 temp[3][3];
	qcd_complex_16 wilsonloopValue;

	qcd_unit3x3(wilsonline);


	for (ll = 0; ll < 2 * l; ll++) //Z 方向
	{
		vp = qcd_LEXIC( (lt - ll - 1 + geo->L[0]) % geo->L[0], lx, ly, lz , geo->lL);
		qcd_mulAdjoint3x3(temp, wilsonline, u->D[vp][0]);
		qcd_copy3x3(wilsonline, temp);

	}
	// 暂时 bprep 只能为正
	for (lbprep = 0; lbprep < bprep; lbprep++)
	{
		vp = qcd_LEXIC((lt - l * 2 + geo->L[0]) % geo->L[0], (lx + lbprep + geo->L[1]) % geo->L[1], ly, lz , geo->lL);
		qcd_mul3x3(temp, wilsonline, u->D[vp][1]);
		qcd_copy3x3(wilsonline, temp);
	}

	for (ll = 0; ll < 2 * l; ll++)
	{
		vp = qcd_LEXIC( (lt - l * 2 + ll + geo->L[0]) % geo->L[0], (lx + bprep + geo->L[1]) % geo->L[1], ly, lz , geo->lL);
		qcd_mul3x3(temp, wilsonline, u->D[vp][0]);
		qcd_copy3x3(wilsonline, temp);
	}


	for (lbprep = 0; lbprep < bprep; lbprep++)
	{
		vp = qcd_LEXIC(lt, (lx + bprep - 1 - lbprep + geo->L[1]) % geo->L[1], ly, lz, geo->lL);
		qcd_mulAdjoint3x3(temp, wilsonline, u->D[vp][1]);
		qcd_copy3x3(wilsonline, temp);
	}

	qcd_dagger3x3(wilsonline);
	wilsonloopValue = qcd_trace3x3(wilsonline);
	return wilsonloopValue;
}*/








void pion_quasi_TMDWF_soft_function(qcd_propagator *dprop, qcd_geometry *geo,qcd_gaugeField *u,int x_src[4], FILE *fp_corr_p)
{
	qcd_int_8 v3, v, vp;
	qcd_int_4 lx, ly, lz,lt;
	qcd_int_4 x, y, z, t;
	qcd_int_4 mu, nu, ku,lu;
	qcd_int_4 c1, c2, c1p, c2p;
	qcd_int_4 ctr = 0;
	qcd_int_4 ctr2;
	qcd_int_2 g5cg0g5_ind[16 * 16][3];
	qcd_complex_16 g5cg0g5_val[16 * 16];
	qcd_complex_16 wilsonline[3][3];
	qcd_complex_16 *block;
	qcd_complex_16 *block2;
	qcd_complex_16 corr, corr2;
	qcd_complex_16 corrw, corrw2; //真空期望值
	qcd_complex_16 C,C2;
	qcd_real_8 tmp;
	qcd_int_4 l, bprep;
	qcd_int_4 mom[3];


// 为了简化 z先取为零

	for (mu = 0; mu < 4; mu++)
		for (ku = 0; ku < 4; ku++)
			for (nu = 0; nu < 4; nu++)
			{
				C = qcd_CMUL(qcd_GAMMA[5][mu][ku],qcd_GAMMAG5[0][ku][nu]);
				//C = qcd_CMUL(qcd_GAMMA[5][mu][ku], qcd_GAMMA[5][ku][nu]);
				if (qcd_NORM(C) > 1e-3)
				{
					g5cg0g5_val[ctr].re = C.re;
					g5cg0g5_val[ctr].im = C.im;
					g5cg0g5_ind[ctr][0] = mu;
					g5cg0g5_ind[ctr][1] = nu;
					g5cg0g5_ind[ctr][2] = ku;
					ctr++;
				}
			}

	block = (qcd_complex_16*)malloc(geo->lV3 * sizeof(qcd_complex_16));
	block2 = (qcd_complex_16*)malloc(geo->lV3 * sizeof(qcd_complex_16));

	for (bprep = -12; bprep <= 12; bprep++)
	{
		for (l = 0; l <= 12; l++)
		{
			for (v3 = 0; v3 < geo->lV3; v3++)   //set blocks to zero
				block[v3] = (qcd_complex_16) { 0, 0 };

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
							create_wilson_line(wilsonline, lt, lx, ly, lz, l, bprep, geo, u);
							block2[v3] = create_wilson_loop(lt, lx, ly, lz, 2*l, abs(bprep), geo, u);

							for (ctr2 = 0; ctr2 < ctr; ctr2++)
							{
								
								for (lu = 0; lu < 4; lu++)
								{
									for (c1 = 0; c1 < 3; c1++)
										for (c1p = 0; c1p < 3; c1p++)
											for (c2p = 0; c2p < 3; c2p++)
											{
												mu = g5cg0g5_ind[ctr2][0];
												nu = g5cg0g5_ind[ctr2][1];
												ku = g5cg0g5_ind[ctr2][2];

												block[v3] = qcd_CADD(block[v3],
													qcd_CMUL(qcd_CMUL(qcd_CONJ(dprop->D[vp][mu][lu][c2p][c1p]), wilsonline[c2p][c1]),
														qcd_CMUL(g5cg0g5_val[ctr2], dprop->D[v][nu][lu][c1][c1p])));

											//block[v3] = qcd_CADD(block[v3],
											//	qcd_CMUL(qcd_CONJ(dprop->D[vp][mu][lu][c1][c1p]), dprop->D[v][mu][lu][c1][c1p]));
														 



											}
								}

							}
						}

			}


			for (mom[2] = 0; mom[2] < 3; mom[2]++)

			{
				mom[0] = 0;
				mom[1] = 0;
				if (geo->myid == 0)
				{
					// fprintf(fp_corr_p, "%+i %+i %+i %+i %+i", wl, mom[j][0], mom[j][1], mom[j][2],wl);
					fprintf(fp_corr_p, " %d  %d  %d  %d  %d ", bprep,l,mom[0],mom[1],mom[2]);

				}
				//if (geo->myid == 0) printf("%d %d %d %d \n", x_src[0], x_src[1], x_src[2], x_src[3]);
				corr = (qcd_complex_16) { 0, 0 };
				corrw = (qcd_complex_16) { 0, 0 };

				if (lt >= 0 && lt < geo->lL[0])  //inside the local lattice, otherwise nothing to calculate
				{
					for (lx = 0; lx < geo->lL[1]; lx++)
						for (ly = 0; ly < geo->lL[2]; ly++)
							for (lz = 0; lz < geo->lL[3]; lz++)
							{
								v3 = qcd_LEXIC0(lx, ly, lz, geo->lL);
								x = lx + geo->Pos[1] * geo->lL[1] - x_src[1];
								y = ly + geo->Pos[2] * geo->lL[2] - x_src[2];
								z = lz + geo->Pos[3] * geo->lL[3] - x_src[3];
								tmp = (((double)mom[0] * x) / geo->L[1] + ((double)mom[1] * y) / geo->L[2] + ((double)mom[2] * z) / geo->L[3]) * 2 * M_PI;
								C2 = (qcd_complex_16) { cos(tmp), -sin(tmp) }; //TABULATE FOR LARGE SPEEDUP!!!
								corr = qcd_CADD(corr, qcd_CMUL(block[v3], C2));
								corrw = qcd_CADD(corrw, block2[v3]);
							}
				}

				MPI_Reduce(&(corr.re), &(corr2.re), 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
				MPI_Reduce(&(corrw.re), &(corrw2.re), 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

				if (geo->myid == 0)
				{
					fprintf(fp_corr_p, "%+e %+e %+e %+e\n", corr2.re, corr2.im,corrw2.re,corrw2.im);
				}

			}

		}

	}

	free(block);
	free(block2);
}






void pion_quasi_vacuum_value( qcd_geometry *geo, qcd_gaugeField *u, int x_src[4], FILE *fp_corr_p)
{
	qcd_int_8 v3, v, vp;
	qcd_int_4 lx, ly, lz, lt;
	qcd_int_4 x, y, z, t;
	qcd_int_4 mu, nu, ku, lu;
	qcd_int_4 c1, c2, c1p, c2p;
	qcd_int_4 ctr = 0;
	qcd_int_4 ctr2;
	qcd_int_2 g5cg0g5_ind[16 * 16][3];
	qcd_complex_16 g5cg0g5_val[16 * 16];
	qcd_complex_16 wilsonline[3][3];
	qcd_complex_16 *block;
	qcd_complex_16 *block2;
	qcd_complex_16 corr, corr2;
	qcd_complex_16 corrw, corrw2; //真空期望值
	qcd_complex_16 C, C2;
	qcd_real_8 tmp;
	qcd_int_4 l, bprep;
	qcd_int_4 mom[3];



	block2 = (qcd_complex_16*)malloc(geo->lV3 * sizeof(qcd_complex_16));

	for (bprep = 0; bprep <= 12; bprep++)
	{
		for (l = 0; l <= 12; l++)
		{
			for (v3 = 0; v3 < geo->lV3; v3++)   //set blocks to zero
			{
				block2[v3] = (qcd_complex_16) { 0, 0 };
			}

		//	lt = ((6 + x_src[0]) % geo->L[0]) - geo->Pos[0] * geo->lL[0]; //这里需要更一般化

			for(lt=0;lt<geo->lL[0];lt++)
				for (lx = 0; lx < geo->lL[1]; lx++)
					for (ly = 0; ly < geo->lL[2]; ly++)
						for (lz = 0; lz < geo->lL[3]; lz++)
						{
							v = qcd_LEXIC(lt, lx, ly, lz, geo->lL);
							vp = qcd_LEXIC(lt, (lx + bprep + geo->L[1]) % geo->L[1], ly, lz, geo->lL);
							v3 = qcd_LEXIC0(lx, ly, lz, geo->lL);
							
							block2[v3] = qcd_CADD(block2[v3],create_wilson_loop(lt, lx, ly, lz, 2*l, abs(bprep), geo, u));
							block2[v3] = qcd_CADD(block2[v3], create_wilson_loop(lt, lx, ly, lz,  abs(bprep), 2 * l, geo, u));
							
							
						}

			
				//if (geo->myid == 0) printf("%d %d %d %d \n", x_src[0], x_src[1], x_src[2], x_src[3]);
				corrw = (qcd_complex_16) { 0, 0 };



				
					for (lx = 0; lx < geo->lL[1]; lx++)
						for (ly = 0; ly < geo->lL[2]; ly++)
							for (lz = 0; lz < geo->lL[3]; lz++)
							{
								v3 = qcd_LEXIC0(lx, ly, lz, geo->lL);
								corrw = qcd_CADD(corrw, block2[v3]);
							}
				

				MPI_Reduce(&(corrw.re), &(corrw2.re), 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
				

				if (geo->myid == 0)
				{
					fprintf(fp_corr_p, "%d %d %+e %+e\n", bprep, l,corrw2.re/geo->V/2/3, corrw2.im/ geo->V/2/3);
				}

			

		}

	}

	free(block2);
}


