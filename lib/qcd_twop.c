#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <qcd.h>
#include <qcd_projectors.h>
#include<time.h>


void twop_pi0_disconnect(qcd_vector * uprop, qcd_vector * dprop, qcd_geometry *geo, qcd_vector *vecsou, qcd_vector *vecsou2, \
	qcd_int_4 x_src[4], qcd_int_2 t_start, qcd_int_2 t_stop, qcd_int_4 mom[3], FILE *fp_corr_p)
{
	qcd_int_4 ctr, ctr2, cc1, cc2, lx, ly, lz, x, y, z, lx2, ly2, lz2;
	qcd_int_4 c1, c2, c3, c1p, c2p, c3p;
	qcd_int_4 mu, nu, ku, lu, gi, al, be;
	qcd_int_4 lx_src[4], i;
	qcd_int_8 v3, v, v2, vtemp;

	qcd_real_8 tmp;
	qcd_complex_16 *block;
	qcd_complex_16 C, C2, corr, corr2, corrdisconn, corrvacuum1, corrvacuum2;
	qcd_complex_16 block2, vacuum1, vacuum2;
	qcd_complex_16 imag = { 0, 1 };
	qcd_complex_16 dpropsou[4][3], upropsou[4][3], vecsousou[4][3], vecsou2sou[4][3];

	qcd_int_4 t, lt;


	block = (qcd_complex_16*)malloc(geo->lV3 * sizeof(qcd_complex_16));

	for (i = 0; i < 4; i++)
		lx_src[i] = x_src[i] - geo->Pos[i] * geo->lL[i];


	if ((lx_src[0] >= 0) && (lx_src[0] < geo->lL[0]) && (lx_src[1] >= 0) && (lx_src[1] < geo->lL[1]) && \
		(lx_src[2] >= 0) && (lx_src[2] < geo->lL[2]) && (lx_src[3] >= 0) && (lx_src[3] < geo->lL[3]))
	{
		printf("source is in myid %i", geo->myid);
		vtemp = qcd_LEXIC(lx_src[0], lx_src[1], lx_src[2], lx_src[3], geo->lL);
		qcd_copy3x3(upropsou, uprop->D[vtemp]);
		qcd_copy3x3(dpropsou, dprop->D[vtemp]);
		qcd_copy3x3(vecsousou, vecsou->D[vtemp]);
		qcd_copy3x3(vecsou2sou, vecsou2->D[vtemp]);
	}

	MPI_Bcast(&(upropsou[0][0].re), 24, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&(dpropsou[0][0].re), 24, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&(vecsousou[0][0].re), 24, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&(vecsou2sou[0][0].re), 24, MPI_DOUBLE, 0, MPI_COMM_WORLD);



	for (t = t_start; t <= t_stop; t++)
	{


		lt = ((t + x_src[0]) % geo->L[0]) - geo->Pos[0] * geo->lL[0];

		if (geo->myid == 0) printf("t=%i\n", t);

		for (v3 = 0; v3 < geo->lV3; v3++)   //set blocks to zero
		{
			block[v3] = (qcd_complex_16) { 0, 0 };

		}
		block2 = (qcd_complex_16) { 0, 0 };
		vacuum1 = (qcd_complex_16) { 0, 0 };
		vacuum2 = (qcd_complex_16) { 0, 0 };

		if (lt >= 0 && lt < geo->lL[0])  //inside the local lattice, otherwise nothing to calculate
		{

			for (be = 0; be < 4; be++)
				for (c1p = 0; c1p < 3; c1p++)
				{
					vacuum2 = qcd_CADD(vacuum2, qcd_CADD(qcd_CMUL(upropsou[be][c1p], qcd_CONJ(vecsousou[be][c1p])),
						qcd_CMUL(dpropsou[be][c1p], qcd_CONJ(vecsou2sou[be][c1p]))));
					for (al = 0; al < 4; al++)
						for (c1 = 0; c1 < 3; c1++)
						{

							for (lx = 0; lx < geo->lL[1]; lx++)
								for (ly = 0; ly < geo->lL[2]; ly++)
									for (lz = 0; lz < geo->lL[3]; lz++)

									{

										v3 = qcd_LEXIC0(lx, ly, lz, geo->lL);
										v = qcd_LEXIC(lt, lx, ly, lz, geo->lL);
										//v2 = qcd_LEXIC(((x_src[0]) % geo->L[0]) - geo->Pos[0] * geo->lL[0], lx2, ly2, lz2, geo->lL)

										block[v3] = qcd_CADD(block[v3]
											, qcd_CMUL(qcd_CMUL(uprop->D[v][al][c1], qcd_CONJ(vecsousou[be][c1p])),
												qcd_CMUL(upropsou[be][c1p], qcd_CONJ(vecsou->D[v][al][c1]))));

										block[v3] = qcd_CADD(block[v3]
											, qcd_CMUL(qcd_CMUL(dprop->D[v][al][c1], qcd_CONJ(vecsou2sou[be][c1p])),
												qcd_CMUL(dpropsou[be][c1p], qcd_CONJ(vecsou2->D[v][al][c1]))));


										block2 = qcd_CADD(block2
											, qcd_CMUL(qcd_CADD(qcd_CMUL(uprop->D[v][al][c1], qcd_CONJ(vecsou->D[v][al][c1])),
												qcd_CMUL(dprop->D[v][al][c1], qcd_CONJ(vecsou2->D[v][al][c1]))),
												qcd_CADD(qcd_CMUL(upropsou[be][c1p], qcd_CONJ(vecsousou[be][c1p])),
													qcd_CMUL(dpropsou[be][c1p], qcd_CONJ(vecsou2sou[be][c1p])))));


										vacuum1 = qcd_CADD(vacuum1, qcd_CADD(qcd_CMUL(uprop->D[v][al][c1], qcd_CONJ(vecsou->D[v][al][c1])),
											qcd_CMUL(dprop->D[v][al][c1], qcd_CONJ(vecsou2->D[v][al][c1]))));





									}
							//nonvanishing projector condition
						}
				}
			//Fourier transform time-slice
		}


		if ((geo->myid) == 0)
		{
			fprintf(fp_corr_p, "%i %+i %+i %+i ", t, mom[0], mom[1], mom[2]);
		}

		corr = (qcd_complex_16) { 0, 0 };

		if (lt >= 0 && lt < geo->lL[0])
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
					}
		}
		//printf("process %i: corr = %f %+fi\n",myid,0.5*corr.re,0.5*corr.im);
		MPI_Reduce(&(corr.re), &(corr2.re), 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&(block2.re), &(corrdisconn.re), 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&(vacuum1.re), &(corrvacuum1.re), 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&(vacuum2.re), &(corrvacuum2.re), 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		if (geo->myid == 0)
		{
			fprintf(fp_corr_p, "%+e %+e %+e %+e %+e %+e %+e %+e \n", corr2.re, corr2.im, corrdisconn.re, corrdisconn.im, corrvacuum1.re, corrvacuum1.im, corrvacuum2.re, corrvacuum2.im);
		}


		//end lt inside local block condition


	}

}


void twop_delta(qcd_propagator * uprop, qcd_propagator* dprop, qcd_geometry *geo, qcd_int_4 x_src[4], qcd_int_2 t_start, qcd_int_2 t_stop, qcd_int_4 mom[3], FILE *fp_corr_p)
{
	qcd_int_4 ctr, ctr2, cc1, cc2, lx, ly, lz, x, y, z;
	qcd_int_4 c1, c2, c3, c1p, c2p, c3p;
	qcd_int_4 mu, nu, ku, lu, gi, al, be;
	qcd_int_8 v3, v;

	qcd_real_8 tmp;
	qcd_complex_16 *block;
	qcd_complex_16 C, C2, corr, corr2;

	qcd_int_4 t, lt;
	qcd_int_2 cgicgi_ind[16 * 16][4];
	qcd_complex_16 cgicgi_val[16 * 16];

	block = (qcd_complex_16*)malloc(geo->lV3 * sizeof(qcd_complex_16));

	ctr = 0;
	for (mu = 0; mu < 4; mu++)
		for (nu = 0; nu < 4; nu++)
			for (ku = 0; ku < 4; ku++)
				for (lu = 0; lu < 4; lu++)
					for (gi = 1; gi < 4; gi++)
					{

						C = qcd_CMUL(qcd_CGAMMA[gi][mu][nu], qcd_BAR_CGAMMA[gi][ku][lu]);
						if (qcd_NORM(C) > 1e-3)
						{
							cgicgi_val[ctr].re = C.re;
							cgicgi_val[ctr].im = C.im;
							cgicgi_ind[ctr][0] = mu;
							cgicgi_ind[ctr][1] = nu;
							cgicgi_ind[ctr][2] = ku;
							cgicgi_ind[ctr][3] = lu;
							ctr++;
						}
					}

	//printf("ID %d start corr", geo->myid);
	for (t = t_start; t <= t_stop; t++)
	{
		lt = ((t + x_src[0]) % geo->L[0]) - geo->Pos[0] * geo->lL[0];

		if (geo->myid == 0) printf("t=%i\n", t);
		for (v3 = 0; v3 < geo->lV3; v3++)   //set blocks to zero
			block[v3] = (qcd_complex_16) { 0, 0 };

		if (lt >= 0 && lt < geo->lL[0])  //inside the local lattice, otherwise nothing to calculate
		{
			for (al = 0; al < 4; al++)
				for (be = 0; be < 4; be++)
					if (qcd_NORM(qcd_ONE_PLUS_GAMMA[0][al][be]) > 0.0001)
					{
						/* alpha-beta-component*/
						//if (myid == 0) printf("process %i: calculating %i-%i component\n", myid, al, be);
						for (ctr2 = 0; ctr2 < ctr; ctr2++)
						{
							mu = cgicgi_ind[ctr2][0];
							nu = cgicgi_ind[ctr2][1];
							ku = cgicgi_ind[ctr2][2];
							lu = cgicgi_ind[ctr2][3];
							for (cc1 = 0; cc1 < 6; cc1++)
							{
								c1 = qcd_EPS[cc1][0];
								c2 = qcd_EPS[cc1][1];
								c3 = qcd_EPS[cc1][2];
								for (cc2 = 0; cc2 < 6; cc2++)
								{
									c1p = qcd_EPS[cc2][0];
									c2p = qcd_EPS[cc2][1];
									c3p = qcd_EPS[cc2][2];
									for (lx = 0; lx < geo->lL[1]; lx++)
										for (ly = 0; ly < geo->lL[2]; ly++)
											for (lz = 0; lz < geo->lL[3]; lz++)
											{
												v3 = qcd_LEXIC0(lx, ly, lz, geo->lL);
												v = qcd_LEXIC(lt, lx, ly, lz, geo->lL);

												block[v3] = qcd_CSUB(block[v3]
													, qcd_CMUL(qcd_ONE_PLUS_GAMMA[0][al][be]
														, qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(cgicgi_val[ctr2]
															, uprop->D[v][mu][ku][c1][c1p])
															, dprop->D[v][nu][lu][c2][c2p])
															, uprop->D[v][be][al][c3][c3p])
															, qcd_SGN_EPS[cc1] * qcd_SGN_EPS[cc2] * 4)));

												block[v3] = qcd_CSUB(block[v3]
													, qcd_CMUL(qcd_ONE_PLUS_GAMMA[0][al][be]
														, qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(cgicgi_val[ctr2]
															, uprop->D[v][mu][al][c1][c1p])
															, dprop->D[v][nu][lu][c2][c2p])
															, uprop->D[v][be][ku][c3][c3p])
															, qcd_SGN_EPS[cc1] * qcd_SGN_EPS[cc2] * 4)));



												block[v3] = qcd_CSUB(block[v3]
													, qcd_CMUL(qcd_ONE_PLUS_GAMMA[0][al][be]
														, qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(cgicgi_val[ctr2]
															, uprop->D[v][mu][ku][c1][c1p])
															, dprop->D[v][nu][al][c2][c3p])
															, uprop->D[v][be][lu][c3][c2p])
															, qcd_SGN_EPS[cc1] * qcd_SGN_EPS[cc2] * (-2))));
												block[v3] = qcd_CSUB(block[v3]
													, qcd_CMUL(qcd_ONE_PLUS_GAMMA[0][al][be]
														, qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(cgicgi_val[ctr2]
															, uprop->D[v][mu][lu][c1][c1p])
															, dprop->D[v][nu][al][c2][c3p])
															, uprop->D[v][be][ku][c3][c2p])
															, qcd_SGN_EPS[cc1] * qcd_SGN_EPS[cc2] * (-2))));

												block[v3] = qcd_CSUB(block[v3]
													, qcd_CMUL(qcd_ONE_PLUS_GAMMA[0][al][be]
														, qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(cgicgi_val[ctr2]
															, uprop->D[v][mu][ku][c1][c1p])
															, dprop->D[v][be][lu][c3][c2p])
															, uprop->D[v][nu][al][c2][c3p])
															, qcd_SGN_EPS[cc1] * qcd_SGN_EPS[cc2] * (-2))));
												block[v3] = qcd_CSUB(block[v3]
													, qcd_CMUL(qcd_ONE_PLUS_GAMMA[0][al][be]
														, qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(cgicgi_val[ctr2]
															, uprop->D[v][mu][al][c1][c1p])
															, dprop->D[v][be][lu][c3][c2p])
															, uprop->D[v][nu][ku][c2][c3p])
															, qcd_SGN_EPS[cc1] * qcd_SGN_EPS[cc2] * (-2))));

												block[v3] = qcd_CSUB(block[v3]
													, qcd_CMUL(qcd_ONE_PLUS_GAMMA[0][al][be]
														, qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(cgicgi_val[ctr2]
															, uprop->D[v][mu][ku][c1][c1p])
															, dprop->D[v][be][al][c3][c3p])
															, uprop->D[v][nu][lu][c2][c2p])
															, qcd_SGN_EPS[cc1] * qcd_SGN_EPS[cc2])));
												block[v3] = qcd_CSUB(block[v3]
													, qcd_CMUL(qcd_ONE_PLUS_GAMMA[0][al][be]
														, qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(cgicgi_val[ctr2]
															, uprop->D[v][mu][lu][c1][c1p])
															, dprop->D[v][be][al][c3][c3p])
															, uprop->D[v][nu][ku][c2][c2p])
															, qcd_SGN_EPS[cc1] * qcd_SGN_EPS[cc2])));

											}//space loop
								}//color2 loop    
							}			//color1 loop
						}//nonvanishing cg5cg5 loop
					}//nonvanishing projector condition
		}
		//Fourier transform time-slice

		//if (geo->myid == 0) printf("now fourier trans\n");
//for (j = 0; j < momi; j++)
//{
		if (geo->myid == 0)
		{
			fprintf(fp_corr_p, "%i %+i %+i %+i ", t, mom[0], mom[1], mom[2]);
		}
		corr = (qcd_complex_16) { 0, 0 };
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
					}
		}
		//if (geo->myid == 0) printf("now fourier trans finish \n");
		//printf("process %i: corr = %f %+fi\n",myid,0.5*corr.re,0.5*corr.im);
		MPI_Reduce(&(corr.re), &(corr2.re), 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		if (geo->myid == 0)
		{
			fprintf(fp_corr_p, "%+e %+e\n", corr2.re*0.25, corr2.im*0.25);   //投影算符的1/4
		}
		//if (geo->myid == 0) printf("one corr finish\n");
	}

	//end lt inside local block condition

	free(block);


}


void twop_proton(qcd_propagator * uprop, qcd_propagator* dprop, qcd_geometry *geo, qcd_int_4 x_src[4], qcd_int_2 t_start, qcd_int_2 t_stop, qcd_int_4 mom[3], FILE *fp_corr_p)
{
	qcd_int_4 ctr, ctr2, cc1, cc2, lx, ly, lz, x, y, z;
	qcd_int_4 c1, c2, c3, c1p, c2p, c3p;
	qcd_int_4 mu, nu, ku, lu, gi, al, be;
	qcd_int_8 v3, v;

	qcd_real_8 tmp;
	qcd_complex_16 *block;
	qcd_complex_16 C, C2, corr, corr2;

	qcd_int_4 t, lt;
	qcd_int_2 cgicgi_ind[16 * 16][4];
	qcd_complex_16 cgicgi_val[16 * 16];

	block = (qcd_complex_16*)malloc(geo->lV3 * sizeof(qcd_complex_16));

	ctr = 0;
	for (mu = 0; mu < 4; mu++)
		for (nu = 0; nu < 4; nu++)
			for (ku = 0; ku < 4; ku++)
				for (lu = 0; lu < 4; lu++)
				{

					C = qcd_CMUL(qcd_CGAMMA[5][mu][nu], qcd_BAR_CGAMMA[5][ku][lu]);
					if (qcd_NORM(C) > 1e-3)
					{
						cgicgi_val[ctr].re = C.re;
						cgicgi_val[ctr].im = C.im;
						cgicgi_ind[ctr][0] = mu;
						cgicgi_ind[ctr][1] = nu;
						cgicgi_ind[ctr][2] = ku;
						cgicgi_ind[ctr][3] = lu;
						ctr++;
					}
				}

	//printf("ID %d start corr", geo->myid);
	for (t = t_start; t <= t_stop; t++)
	{
		lt = ((t + x_src[0]) % geo->L[0]) - geo->Pos[0] * geo->lL[0];

		if (geo->myid == 0) printf("t=%i\n", t);
		for (v3 = 0; v3 < geo->lV3; v3++)   //set blocks to zero
			block[v3] = (qcd_complex_16) { 0, 0 };

		if (lt >= 0 && lt < geo->lL[0])  //inside the local lattice, otherwise nothing to calculate
		{
			for (al = 0; al < 4; al++)
				for (be = 0; be < 4; be++)
					if (qcd_NORM(qcd_ONE_PLUS_GAMMA[0][al][be]) > 0.0001)
					{
						/* alpha-beta-component*/
						//if (myid == 0) printf("process %i: calculating %i-%i component\n", myid, al, be);
						for (ctr2 = 0; ctr2 < ctr; ctr2++)
						{
							mu = cgicgi_ind[ctr2][0];
							nu = cgicgi_ind[ctr2][1];
							ku = cgicgi_ind[ctr2][2];
							lu = cgicgi_ind[ctr2][3];
							for (cc1 = 0; cc1 < 6; cc1++)
							{
								c1 = qcd_EPS[cc1][0];
								c2 = qcd_EPS[cc1][1];
								c3 = qcd_EPS[cc1][2];
								for (cc2 = 0; cc2 < 6; cc2++)
								{
									c1p = qcd_EPS[cc2][0];
									c2p = qcd_EPS[cc2][1];
									c3p = qcd_EPS[cc2][2];
									for (lx = 0; lx < geo->lL[1]; lx++)
										for (ly = 0; ly < geo->lL[2]; ly++)
											for (lz = 0; lz < geo->lL[3]; lz++)
											{
												v3 = qcd_LEXIC0(lx, ly, lz, geo->lL);
												v = qcd_LEXIC(lt, lx, ly, lz, geo->lL);

												block[v3] = qcd_CSUB(block[v3]
													, qcd_CMUL(qcd_ONE_PLUS_GAMMA[0][al][be]
														, qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(cgicgi_val[ctr2]
															, uprop->D[v][mu][ku][c1][c1p])
															, dprop->D[v][nu][lu][c2][c2p])
															, uprop->D[v][be][al][c3][c3p])
														, qcd_SGN_EPS[cc1] * qcd_SGN_EPS[cc2])));

												block[v3] = qcd_CSUB(block[v3]
													, qcd_CMUL(qcd_ONE_PLUS_GAMMA[0][al][be]
														,qcd_CSCALE( qcd_CMUL(qcd_CMUL(qcd_CMUL(cgicgi_val[ctr2]
															, uprop->D[v][mu][al][c1][c1p])
															, dprop->D[v][nu][lu][c2][c2p])
															, uprop->D[v][be][ku][c3][c3p])
														, qcd_SGN_EPS[cc1] * qcd_SGN_EPS[cc2])));



											}//space loop
								}//color2 loop    
							}			//color1 loop
						}//nonvanishing cg5cg5 loop
					}//nonvanishing projector condition
		}
		//Fourier transform time-slice

		//if (geo->myid == 0) printf("now fourier trans\n");
//for (j = 0; j < momi; j++)
//{
		if (geo->myid == 0)
		{
			fprintf(fp_corr_p, "%i %+i %+i %+i ", t, mom[0], mom[1], mom[2]);
		}
		corr = (qcd_complex_16) { 0, 0 };
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
					}
		}
		//if (geo->myid == 0) printf("now fourier trans finish \n");
		//printf("process %i: corr = %f %+fi\n",myid,0.5*corr.re,0.5*corr.im);
		MPI_Reduce(&(corr.re), &(corr2.re), 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		if (geo->myid == 0)
		{
			fprintf(fp_corr_p, "%+e %+e\n", corr2.re*0.25, corr2.im*0.25);  //投影算符中漏的1/4
		}
		//if (geo->myid == 0) printf("one corr finish\n");
	}

	//end lt inside local block condition

	free(block);


}