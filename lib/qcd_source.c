
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <qcd.h>
#include <qcd_projectors.h>
#include<time.h>






void mksource_su_proton(qcd_vector *souvec, qcd_propagator *uprop, qcd_propagator *dprop, qcd_gaugeField *uAPE, qcd_geometry *geo, \
	qcd_int_4 x_src[4], qcd_int_4 is, qcd_real_8 alpha, qcd_uint_4 nsmear, qcd_real_8 keci, qcd_int_4 mom[3], qcd_int_4 tsink)
{
	qcd_int_4 mu, nu, ku, lu, be, al, ta;
	qcd_int_4 c1, c2, c3, c1p, c2p, c3p, cc1, cc2, a, b;
	qcd_int_4 ctr, ctr2, row, col, i, gi, v;
	qcd_int_4 x, y, z, lx, ly, lz, t, lt;
	qcd_int_4 numb_sources = 12;

	qcd_real_8 tmp;
	qcd_complex_16 C, C2;
	qcd_complex_16 imag = { 0, 1 };

	qcd_int_2 cgicgi_ind[16 * 16 * 16][6];
	qcd_complex_16 cgicgi_val[16 * 16 * 16];


	qcd_propagator soupro;

	numb_sources = 12;


	ctr = 0;
	for (mu = 0; mu < 4; mu++)
		for (nu = 0; nu < 4; nu++)
			for (ku = 0; ku < 4; ku++)
				for (lu = 0; lu < 4; lu++)
					for (be = 0; be < 4; be++)
						for (al = 0; al < 4; al++)
						{
							C.re = 0;
							C.im = 0;
							C = qcd_CADD(C, qcd_CMUL(PROJECTOR[17][be][al], qcd_CMUL(qcd_CGAMMA[5][mu][nu], qcd_BAR_CGAMMA[5][ku][lu])));

							C = qcd_CADD(C, qcd_CMUL(PROJECTOR[17][lu][al], qcd_CMUL(qcd_CGAMMA[5][mu][nu], qcd_BAR_CGAMMA[5][ku][be])));

							C = qcd_CADD(C, qcd_CMUL(PROJECTOR[17][be][mu], qcd_CMUL(qcd_CGAMMA[5][al][nu], qcd_BAR_CGAMMA[5][ku][lu])));

							C = qcd_CADD(C, qcd_CMUL(PROJECTOR[17][lu][mu], qcd_CMUL(qcd_CGAMMA[5][al][nu], qcd_BAR_CGAMMA[5][ku][be])));


							if (qcd_NORM(C) > 1e-3)
							{
								cgicgi_val[ctr].re = C.re;
								cgicgi_val[ctr].im = C.im;
								cgicgi_ind[ctr][0] = mu;
								cgicgi_ind[ctr][1] = nu;
								cgicgi_ind[ctr][2] = ku;
								cgicgi_ind[ctr][3] = lu;
								cgicgi_ind[ctr][4] = al;
								cgicgi_ind[ctr][5] = be;
								ctr++;
							}
						}
	qcd_initPropagator(&soupro, geo);
	qcd_zeroPropagator(&soupro);
	//qcd_zeroVector(souvec);
	t = tsink;
	lt = ((t + x_src[0]) % geo->L[0]) - geo->Pos[0] * geo->lL[0];
	if (lt >= 0 && lt < geo->lL[0])  //inside the local lattice, otherwise nothing to calculate
	{

		for (ctr2 = 0; ctr2 < ctr; ctr2++)
		{
			mu = cgicgi_ind[ctr2][0];
			nu = cgicgi_ind[ctr2][1];
			ku = cgicgi_ind[ctr2][2];
			lu = cgicgi_ind[ctr2][3];
			al = cgicgi_ind[ctr2][4];
			be = cgicgi_ind[ctr2][5];
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
					for (ta = 0; ta < 4; ta++)
					{
						if (qcd_NORM(qcd_CADD(qcd_GAMMA[5][ta][mu], qcd_CMUL(imag, PROJECTOR[0][ta][mu]))) > 1e-3)
						{
							for (lx = 0; lx < geo->lL[1]; lx++)
								for (ly = 0; ly < geo->lL[2]; ly++)
									for (lz = 0; lz < geo->lL[3]; lz++)
									{
										//v3 = qcd_LEXIC0(lx, ly, lz, geo.lL);
										v = qcd_LEXIC(lt, lx, ly, lz, geo->lL);

										soupro.D[v][ta][lu][c1][c3p] = qcd_CADD(soupro.D[v][ta][lu][c1][c3p], qcd_CONJ(
											qcd_CMUL(qcd_CADD(qcd_GAMMA[5][ta][mu], qcd_CMUL(imag, PROJECTOR[0][ta][mu])),
												qcd_CSCALE(qcd_CMUL(qcd_CMUL(cgicgi_val[ctr2]
													, dprop->D[v][nu][ku][c2][c2p])
													, uprop->D[v][al][be][c3][c1p])
													, qcd_SGN_EPS[cc1] * qcd_SGN_EPS[cc2]))));
									}
						}
					}//space loop
				}//color2 loop    
			}//color1 loop
		}//nonvanishing cg5cg5 loop
	//nonvanishing projector condition

	//Fourier transform time-slice


		for (lx = 0; lx < geo->lL[1]; lx++)
			for (ly = 0; ly < geo->lL[2]; ly++)
				for (lz = 0; lz < geo->lL[3]; lz++)
				{
					v = qcd_LEXIC(lt, lx, ly, lz, geo->lL);
					x = lx + geo->Pos[1] * geo->lL[1] - x_src[1];
					y = ly + geo->Pos[2] * geo->lL[2] - x_src[2];
					z = lz + geo->Pos[3] * geo->lL[3] - x_src[3];
					tmp = (((double)mom[0] * x) / geo->L[1] + ((double)mom[1] * y) / geo->L[2] + ((double)mom[2] * z) / geo->L[3]) * 2 * M_PI;
					C2 = (qcd_complex_16) { cos(tmp), sin(tmp) }; //TABULATE FOR LARGE SPEEDUP!!!
					for (mu = 0; mu < 4; mu++)
						for (nu = 0; nu < 4; nu++)
							for (a = 0; a < 3; a++)
								for (b = 0; b < 3; b++)
									soupro.D[v][mu][nu][a][b] = qcd_CMUL(soupro.D[v][mu][nu][a][b], C2);
				}
	}


	row = is / 3;
	col = is % 3;
	qcd_copyVectorPropagator(souvec, &soupro, row, col);


	for (i = 0; i < nsmear; i++)
	{
		if (qcd_momgaussIteration3d(souvec, uAPE, alpha, (t + x_src[0]) % geo->L[0], mom, keci))
		{
			fprintf(stderr, "process %i: Error while smearing!\n", geo->myid);
			exit(EXIT_FAILURE);
		}
	}

	if (geo->myid == 0) printf(" Done source_su vector: %4d / %4d\n", is, numb_sources);


	qcd_destroyPropagator(&soupro);
}







void mksource_sd_proton(qcd_vector *souvec, qcd_propagator *uprop, qcd_propagator *dprop, qcd_gaugeField *uAPE, qcd_geometry *geo, \
	qcd_int_4 x_src[4], qcd_int_4 is, qcd_real_8 alpha, qcd_uint_4 nsmear, qcd_real_8 keci, qcd_int_4 mom[3], qcd_int_4 tsink)
{

	qcd_int_4 mu, nu, ku, lu, be, al, ta;
	qcd_int_4 c1, c2, c3, c1p, c2p, c3p, cc1, cc2, a, b;
	qcd_int_4 ctr, ctr2, row, col, i, gi, v;
	qcd_int_4 x, y, z, lx, ly, lz, t, lt;
	qcd_int_4 numb_sources = 12;

	qcd_real_8 tmp;
	qcd_complex_16 C, C2;
	qcd_complex_16 imag = { 0, 1 };

	qcd_int_2 cgicgi_ind[16 * 16 * 16][6];
	qcd_complex_16 cgicgi_val[16 * 16 * 16];


	qcd_propagator soupro;



	numb_sources = 12;

	ctr = 0;
	for (mu = 0; mu < 4; mu++)
		for (nu = 0; nu < 4; nu++)
			for (ku = 0; ku < 4; ku++)
				for (lu = 0; lu < 4; lu++)
					for (be = 0; be < 4; be++)
						for (al = 0; al < 4; al++)
						{
							C.re = 0;
							C.im = 0;
							C = qcd_CADD(C, qcd_CMUL(PROJECTOR[17][be][al], qcd_CMUL(qcd_CGAMMA[5][mu][nu], qcd_BAR_CGAMMA[5][ku][lu])));
							if (qcd_NORM(C) > 1e-3)
							{
								cgicgi_val[ctr].re = C.re;
								cgicgi_val[ctr].im = C.im;
								cgicgi_ind[ctr][0] = mu;
								cgicgi_ind[ctr][1] = nu;
								cgicgi_ind[ctr][2] = ku;
								cgicgi_ind[ctr][3] = lu;
								cgicgi_ind[ctr][4] = al;
								cgicgi_ind[ctr][5] = be;
								ctr++;
							}
						}
	qcd_initPropagator(&soupro, geo);
	qcd_zeroPropagator(&soupro);
	//qcd_zeroVector(souvec);
	t = tsink;
	lt = ((t + x_src[0]) % geo->L[0]) - geo->Pos[0] * geo->lL[0];
	if (lt >= 0 && lt < geo->lL[0])  //inside the local lattice, otherwise nothing to calculate
	{
		if (geo->myid == 0) printf("t=%i\n", t);

		//if (myid == 0) printf("process %i: calculating %i-%i component\n", myid, al, be);
		for (ctr2 = 0; ctr2 < ctr; ctr2++)
		{
			mu = cgicgi_ind[ctr2][0];
			nu = cgicgi_ind[ctr2][1];
			ku = cgicgi_ind[ctr2][2];
			lu = cgicgi_ind[ctr2][3];
			al = cgicgi_ind[ctr2][4];
			be = cgicgi_ind[ctr2][5];
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
					for (ta = 0; ta < 4; ta++)
					{
						if (qcd_NORM(qcd_CSUB(qcd_GAMMA[5][ta][nu], qcd_CMUL(imag, PROJECTOR[0][ta][nu]))) > 1e-3)
						{
							for (lx = 0; lx < geo->lL[1]; lx++)
								for (ly = 0; ly < geo->lL[2]; ly++)
									for (lz = 0; lz < geo->lL[3]; lz++)
									{
										//v3 = qcd_LEXIC0(lx, ly, lz, geo.lL);
										v = qcd_LEXIC(lt, lx, ly, lz, geo->lL);

										soupro.D[v][ta][ku][c2][c2p] = qcd_CADD(soupro.D[v][ta][ku][c2][c2p], qcd_CONJ(
											qcd_CMUL(qcd_CSUB(qcd_GAMMA[5][ta][nu], qcd_CMUL(imag, PROJECTOR[0][ta][nu])), qcd_CSCALE(qcd_CMUL(cgicgi_val[ctr2],
												qcd_CSUB(qcd_CMUL(uprop->D[v][mu][lu][c1][c3p], uprop->D[v][al][be][c3][c1p]), qcd_CMUL(uprop->D[v][mu][be][c1][c1p], uprop->D[v][al][lu][c3][c3p])))
												, qcd_SGN_EPS[cc1] * qcd_SGN_EPS[cc2]))));

									}
						}
					}//space loop
				}//color2 loop    
			}//color1 loop
		}//nonvanishing cg5cg5 loop
	//nonvanishing projector condition



		for (lx = 0; lx < geo->lL[1]; lx++)
			for (ly = 0; ly < geo->lL[2]; ly++)
				for (lz = 0; lz < geo->lL[3]; lz++)
				{
					v = qcd_LEXIC(lt, lx, ly, lz, geo->lL);
					x = lx + geo->Pos[1] * geo->lL[1] - x_src[1];
					y = ly + geo->Pos[2] * geo->lL[2] - x_src[2];
					z = lz + geo->Pos[3] * geo->lL[3] - x_src[3];
					tmp = (((double)mom[0] * x) / geo->L[1] + ((double)mom[1] * y) / geo->L[2] + ((double)mom[2] * z) / geo->L[3]) * 2 * M_PI;
					C2 = (qcd_complex_16) { cos(tmp), sin(tmp) }; //TABULATE FOR LARGE SPEEDUP!!!
					for (mu = 0; mu < 4; mu++)
						for (nu = 0; nu < 4; nu++)
							for (a = 0; a < 3; a++)
								for (b = 0; b < 3; b++)
									soupro.D[v][mu][nu][a][b] = qcd_CMUL(soupro.D[v][mu][nu][a][b], C2);

				}

	}

	row = is / 3;
	col = is % 3;

	qcd_copyVectorPropagator(souvec, &soupro, row, col);


	for (i = 0; i < nsmear; i++)
	{
		if (qcd_momgaussIteration3d(souvec, uAPE, alpha, (t + x_src[0]) % geo->L[0], mom, keci))
		{
			fprintf(stderr, "process %i: Error while smearing!\n", geo->myid);
			exit(EXIT_FAILURE);
		}
	}

	if (geo->myid == 0) printf(" Done source_dn vector: %4d / %4d\n", is, numb_sources);

	qcd_destroyPropagator(&soupro);

}







void mksource_su_delta(qcd_vector *souvec, qcd_propagator *uprop, qcd_propagator *dprop, qcd_gaugeField *uAPE, qcd_geometry *geo, \
	qcd_int_4 x_src[4], qcd_int_4 is, qcd_real_8 alpha, qcd_uint_4 nsmear, qcd_real_8 keci, qcd_int_4 mom[3], qcd_int_4 tsink)
{
	qcd_int_4 mu, nu, ku, lu, be, al, ta;
	qcd_int_4 c1, c2, c3, c1p, c2p, c3p, cc1, cc2, a, b;
	qcd_int_4 ctr, ctr2, row, col, i, gi, v;
	qcd_int_4 x, y, z, lx, ly, lz, t, lt;
	qcd_int_4 numb_sources = 12;

	qcd_real_8 tmp;
	qcd_complex_16 C, C2;
	qcd_complex_16 imag = { 0, 1 };

	qcd_int_2 cgicgi_ind[16 * 16 * 16][6];
	qcd_complex_16 cgicgi_val[16 * 16 * 16];


	qcd_propagator soupro;

	numb_sources = 12;


	ctr = 0;
	for (mu = 0; mu < 4; mu++)
		for (nu = 0; nu < 4; nu++)
			for (ku = 0; ku < 4; ku++)
				for (lu = 0; lu < 4; lu++)
					for (be = 0; be < 4; be++)
						for (al = 0; al < 4; al++)
							for (gi = 1; gi < 4; gi++)
							{
								C.re = 0;
								C.im = 0;
								C = qcd_CADD(C, qcd_CSCALE(qcd_CMUL(PROJECTOR[18][be][al], qcd_CMUL(qcd_CGAMMA[gi][mu][nu], qcd_BAR_CGAMMA[gi][ku][lu])), 4));
								C = qcd_CADD(C, qcd_CSCALE(qcd_CMUL(PROJECTOR[18][ku][al], qcd_CMUL(qcd_CGAMMA[gi][mu][nu], qcd_BAR_CGAMMA[gi][be][lu])), 2));
								C = qcd_CADD(C, qcd_CSCALE(qcd_CMUL(PROJECTOR[18][be][nu], qcd_CMUL(qcd_CGAMMA[gi][mu][al], qcd_BAR_CGAMMA[gi][ku][lu])), 2));
								C = qcd_CADD(C, qcd_CSCALE(qcd_CMUL(PROJECTOR[18][ku][nu], qcd_CMUL(qcd_CGAMMA[gi][mu][al], qcd_BAR_CGAMMA[gi][be][lu])), 1));

								C = qcd_CADD(C, qcd_CSCALE(qcd_CMUL(PROJECTOR[18][lu][al], qcd_CMUL(qcd_CGAMMA[gi][mu][nu], qcd_BAR_CGAMMA[gi][ku][be])), 4));
								C = qcd_CADD(C, qcd_CSCALE(qcd_CMUL(PROJECTOR[18][ku][al], qcd_CMUL(qcd_CGAMMA[gi][mu][nu], qcd_BAR_CGAMMA[gi][lu][be])), 2));
								C = qcd_CADD(C, qcd_CSCALE(qcd_CMUL(PROJECTOR[18][lu][nu], qcd_CMUL(qcd_CGAMMA[gi][mu][al], qcd_BAR_CGAMMA[gi][ku][be])), 2));
								C = qcd_CADD(C, qcd_CSCALE(qcd_CMUL(PROJECTOR[18][ku][nu], qcd_CMUL(qcd_CGAMMA[gi][mu][al], qcd_BAR_CGAMMA[gi][lu][be])), 1));

								C = qcd_CADD(C, qcd_CSCALE(qcd_CMUL(PROJECTOR[18][be][mu], qcd_CMUL(qcd_CGAMMA[gi][al][nu], qcd_BAR_CGAMMA[gi][ku][lu])), 4));
								C = qcd_CADD(C, qcd_CSCALE(qcd_CMUL(PROJECTOR[18][ku][mu], qcd_CMUL(qcd_CGAMMA[gi][al][nu], qcd_BAR_CGAMMA[gi][be][lu])), 2));
								C = qcd_CADD(C, qcd_CSCALE(qcd_CMUL(PROJECTOR[18][be][nu], qcd_CMUL(qcd_CGAMMA[gi][al][mu], qcd_BAR_CGAMMA[gi][ku][lu])), 2));
								C = qcd_CADD(C, qcd_CSCALE(qcd_CMUL(PROJECTOR[18][ku][nu], qcd_CMUL(qcd_CGAMMA[gi][al][mu], qcd_BAR_CGAMMA[gi][be][lu])), 1));

								C = qcd_CADD(C, qcd_CSCALE(qcd_CMUL(PROJECTOR[18][lu][mu], qcd_CMUL(qcd_CGAMMA[gi][al][nu], qcd_BAR_CGAMMA[gi][ku][be])), 4));
								C = qcd_CADD(C, qcd_CSCALE(qcd_CMUL(PROJECTOR[18][ku][mu], qcd_CMUL(qcd_CGAMMA[gi][al][nu], qcd_BAR_CGAMMA[gi][lu][be])), 2));
								C = qcd_CADD(C, qcd_CSCALE(qcd_CMUL(PROJECTOR[18][lu][nu], qcd_CMUL(qcd_CGAMMA[gi][al][mu], qcd_BAR_CGAMMA[gi][ku][be])), 2));
								C = qcd_CADD(C, qcd_CSCALE(qcd_CMUL(PROJECTOR[18][ku][nu], qcd_CMUL(qcd_CGAMMA[gi][al][mu], qcd_BAR_CGAMMA[gi][lu][be])), 1));


								if (qcd_NORM(C) > 1e-3)
								{
									cgicgi_val[ctr].re = C.re;
									cgicgi_val[ctr].im = C.im;
									cgicgi_ind[ctr][0] = mu;
									cgicgi_ind[ctr][1] = nu;
									cgicgi_ind[ctr][2] = ku;
									cgicgi_ind[ctr][3] = lu;
									cgicgi_ind[ctr][4] = al;
									cgicgi_ind[ctr][5] = be;
									ctr++;
								}
							}
	qcd_initPropagator(&soupro, geo);
	qcd_zeroPropagator(&soupro);
	//qcd_zeroVector(souvec);
	t = tsink;
	lt = ((t + x_src[0]) % geo->L[0]) - geo->Pos[0] * geo->lL[0];
	if (lt >= 0 && lt < geo->lL[0])  //inside the local lattice, otherwise nothing to calculate
	{

		for (ctr2 = 0; ctr2 < ctr; ctr2++)
		{
			mu = cgicgi_ind[ctr2][0];
			nu = cgicgi_ind[ctr2][1];
			ku = cgicgi_ind[ctr2][2];
			lu = cgicgi_ind[ctr2][3];
			al = cgicgi_ind[ctr2][4];
			be = cgicgi_ind[ctr2][5];
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
					for (ta = 0; ta < 4; ta++)
					{
						if (qcd_NORM(qcd_CADD(qcd_GAMMA[5][ta][mu], qcd_CMUL(imag, PROJECTOR[0][ta][mu]))) > 1e-3)
						{
							for (lx = 0; lx < geo->lL[1]; lx++)
								for (ly = 0; ly < geo->lL[2]; ly++)
									for (lz = 0; lz < geo->lL[3]; lz++)
									{
										//v3 = qcd_LEXIC0(lx, ly, lz, geo.lL);
										v = qcd_LEXIC(lt, lx, ly, lz, geo->lL);

										soupro.D[v][ta][lu][c1][c3p] = qcd_CADD(soupro.D[v][ta][lu][c1][c3p], qcd_CONJ(
											qcd_CMUL(qcd_CADD(qcd_GAMMA[5][ta][mu], qcd_CMUL(imag, PROJECTOR[0][ta][mu])),
												qcd_CSCALE(qcd_CMUL(qcd_CMUL(cgicgi_val[ctr2]
													, dprop->D[v][nu][ku][c2][c2p])
													, uprop->D[v][al][be][c3][c1p])
													, qcd_SGN_EPS[cc1] * qcd_SGN_EPS[cc2]))));
									}
						}
					}//space loop
				}//color2 loop    
			}//color1 loop
		}//nonvanishing cg5cg5 loop
	//nonvanishing projector condition

	//Fourier transform time-slice


		for (lx = 0; lx < geo->lL[1]; lx++)
			for (ly = 0; ly < geo->lL[2]; ly++)
				for (lz = 0; lz < geo->lL[3]; lz++)
				{
					v = qcd_LEXIC(lt, lx, ly, lz, geo->lL);
					x = lx + geo->Pos[1] * geo->lL[1] - x_src[1];
					y = ly + geo->Pos[2] * geo->lL[2] - x_src[2];
					z = lz + geo->Pos[3] * geo->lL[3] - x_src[3];
					tmp = (((double)mom[0] * x) / geo->L[1] + ((double)mom[1] * y) / geo->L[2] + ((double)mom[2] * z) / geo->L[3]) * 2 * M_PI;
					C2 = (qcd_complex_16) { cos(tmp), sin(tmp) }; //TABULATE FOR LARGE SPEEDUP!!!
					for (mu = 0; mu < 4; mu++)
						for (nu = 0; nu < 4; nu++)
							for (a = 0; a < 3; a++)
								for (b = 0; b < 3; b++)
									soupro.D[v][mu][nu][a][b] = qcd_CMUL(soupro.D[v][mu][nu][a][b], C2);
				}
	}


	row = is / 3;
	col = is % 3;
	qcd_copyVectorPropagator(souvec, &soupro, row, col);


	for (i = 0; i < nsmear; i++)
	{
		if (qcd_momgaussIteration3d(souvec, uAPE, alpha, (t + x_src[0]) % geo->L[0], mom, keci))
		{
			fprintf(stderr, "process %i: Error while smearing!\n", geo->myid);
			exit(EXIT_FAILURE);
		}
	}

	if (geo->myid == 0) printf(" Done source_su vector: %4d / %4d\n", is, numb_sources);


	qcd_destroyPropagator(&soupro);
}







void mksource_sd_delta(qcd_vector *souvec, qcd_propagator *uprop, qcd_propagator *dprop, qcd_gaugeField *uAPE, qcd_geometry *geo, \
	qcd_int_4 x_src[4], qcd_int_4 is, qcd_real_8 alpha, qcd_uint_4 nsmear, qcd_real_8 keci, qcd_int_4 mom[3], qcd_int_4 tsink)
{

	qcd_int_4 mu, nu, ku, lu, be, al, ta;
	qcd_int_4 c1, c2, c3, c1p, c2p, c3p, cc1, cc2, a, b;
	qcd_int_4 ctr, ctr2, row, col, i, gi, v;
	qcd_int_4 x, y, z, lx, ly, lz, t, lt;
	qcd_int_4 numb_sources = 12;

	qcd_real_8 tmp;
	qcd_complex_16 C, C2;
	qcd_complex_16 imag = { 0, 1 };

	qcd_int_2 cgicgi_ind[16 * 16 * 16][6];
	qcd_complex_16 cgicgi_val[16 * 16 * 16];


	qcd_propagator soupro;



	numb_sources = 12;

	ctr = 0;
	for (mu = 0; mu < 4; mu++)
		for (nu = 0; nu < 4; nu++)
			for (ku = 0; ku < 4; ku++)
				for (lu = 0; lu < 4; lu++)
					for (be = 0; be < 4; be++)
						for (al = 0; al < 4; al++)
							for (gi = 1; gi < 4; gi++)
							{
								C.re = 0;
								C.im = 0;
								C = qcd_CADD(C, qcd_CSCALE(qcd_CMUL(PROJECTOR[18][be][al], qcd_CMUL(qcd_CGAMMA[gi][mu][nu], qcd_BAR_CGAMMA[gi][ku][lu])), 4));
								C = qcd_CADD(C, qcd_CSCALE(qcd_CMUL(PROJECTOR[18][ku][al], qcd_CMUL(qcd_CGAMMA[gi][mu][nu], qcd_BAR_CGAMMA[gi][be][lu])), 2));
								C = qcd_CADD(C, qcd_CSCALE(qcd_CMUL(PROJECTOR[18][be][nu], qcd_CMUL(qcd_CGAMMA[gi][mu][al], qcd_BAR_CGAMMA[gi][ku][lu])), 2));
								C = qcd_CADD(C, qcd_CSCALE(qcd_CMUL(PROJECTOR[18][ku][nu], qcd_CMUL(qcd_CGAMMA[gi][mu][al], qcd_BAR_CGAMMA[gi][be][lu])), 1));

								if (qcd_NORM(C) > 1e-3)
								{
									cgicgi_val[ctr].re = C.re;
									cgicgi_val[ctr].im = C.im;
									cgicgi_ind[ctr][0] = mu;
									cgicgi_ind[ctr][1] = nu;
									cgicgi_ind[ctr][2] = ku;
									cgicgi_ind[ctr][3] = lu;
									cgicgi_ind[ctr][4] = al;
									cgicgi_ind[ctr][5] = be;
									ctr++;
								}
							}
	qcd_initPropagator(&soupro, geo);
	qcd_zeroPropagator(&soupro);
	//qcd_zeroVector(souvec);
	t = tsink;
	lt = ((t + x_src[0]) % geo->L[0]) - geo->Pos[0] * geo->lL[0];
	if (lt >= 0 && lt < geo->lL[0])  //inside the local lattice, otherwise nothing to calculate
	{
		if (geo->myid == 0) printf("t=%i\n", t);

		//if (myid == 0) printf("process %i: calculating %i-%i component\n", myid, al, be);
		for (ctr2 = 0; ctr2 < ctr; ctr2++)
		{
			mu = cgicgi_ind[ctr2][0];
			nu = cgicgi_ind[ctr2][1];
			ku = cgicgi_ind[ctr2][2];
			lu = cgicgi_ind[ctr2][3];
			al = cgicgi_ind[ctr2][4];
			be = cgicgi_ind[ctr2][5];
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
					for (ta = 0; ta < 4; ta++)
					{
						if (qcd_NORM(qcd_CSUB(qcd_GAMMA[5][ta][nu], qcd_CMUL(imag, PROJECTOR[0][ta][nu]))) > 1e-3)
						{
							for (lx = 0; lx < geo->lL[1]; lx++)
								for (ly = 0; ly < geo->lL[2]; ly++)
									for (lz = 0; lz < geo->lL[3]; lz++)
									{
										//v3 = qcd_LEXIC0(lx, ly, lz, geo.lL);
										v = qcd_LEXIC(lt, lx, ly, lz, geo->lL);

										soupro.D[v][ta][ku][c2][c2p] = qcd_CADD(soupro.D[v][ta][ku][c2][c2p], qcd_CONJ(
											qcd_CMUL(qcd_CSUB(qcd_GAMMA[5][ta][nu], qcd_CMUL(imag, PROJECTOR[0][ta][nu])), qcd_CSCALE(qcd_CMUL(cgicgi_val[ctr2],
												qcd_CSUB(qcd_CMUL(uprop->D[v][mu][lu][c1][c3p], uprop->D[v][al][be][c3][c1p]), qcd_CMUL(uprop->D[v][mu][be][c1][c1p], uprop->D[v][al][lu][c3][c3p])))
												, qcd_SGN_EPS[cc1] * qcd_SGN_EPS[cc2]))));

									}
						}
					}//space loop
				}//color2 loop    
			}//color1 loop
		}//nonvanishing cg5cg5 loop
	//nonvanishing projector condition



		for (lx = 0; lx < geo->lL[1]; lx++)
		for (ly = 0; ly < geo->lL[2]; ly++)
		for (lz = 0; lz < geo->lL[3]; lz++)
		{
			v = qcd_LEXIC(lt, lx, ly, lz, geo->lL);
			x = lx + geo->Pos[1] * geo->lL[1] - x_src[1];
			y = ly + geo->Pos[2] * geo->lL[2] - x_src[2];
			z = lz + geo->Pos[3] * geo->lL[3] - x_src[3];
			tmp = (((double)mom[0] * x) / geo->L[1] + ((double)mom[1] * y) / geo->L[2] + ((double)mom[2] * z) / geo->L[3]) * 2 * M_PI;
			C2 = (qcd_complex_16) { cos(tmp), sin(tmp) }; //TABULATE FOR LARGE SPEEDUP!!!
			for (mu = 0; mu < 4; mu++)
			for (nu = 0; nu < 4; nu++)
			for (a = 0; a < 3; a++)
			for (b = 0; b < 3; b++)
				soupro.D[v][mu][nu][a][b] = qcd_CMUL(soupro.D[v][mu][nu][a][b], C2);

		}

	}

	row = is / 3;
	col = is % 3;

	qcd_copyVectorPropagator(souvec, &soupro, row, col);


	for (i = 0; i < nsmear; i++)
	{
		if (qcd_momgaussIteration3d(souvec, uAPE, alpha, (t + x_src[0]) % geo->L[0], mom, keci))
		{
			fprintf(stderr, "process %i: Error while smearing!\n", geo->myid);
			exit(EXIT_FAILURE);
		}
	}

	if (geo->myid == 0) printf(" Done source_dn vector: %4d / %4d\n", is, numb_sources);

	qcd_destroyPropagator(&soupro);

}








void mksource(qcd_vector *vec, qcd_gaugeField *uAPE, qcd_geometry *geo, qcd_int_4 is, qcd_int_4 x_src[4], qcd_real_8 alpha, qcd_uint_4 nsmear, qcd_real_8 keci, qcd_int_4 mom[3])
{
	qcd_int_4 lx_src[4];
	int numb_sources = 12;
	int i;

	for (i = 0; i < 4; i++)
		lx_src[i] = x_src[i] - geo->Pos[i] * geo->lL[i];

	int mu = is / 3;
	int col = is % 3;
	qcd_zeroVector(vec);
	if ((lx_src[0] >= 0) && (lx_src[0] < geo->lL[0]) && (lx_src[1] >= 0) && (lx_src[1] < geo->lL[1]) && (lx_src[2] >= 0) && (lx_src[2] < geo->lL[2]) && (lx_src[3] >= 0) && (lx_src[3] < geo->lL[3]))
		vec->D[qcd_LEXIC(lx_src[0], lx_src[1], lx_src[2], lx_src[3], geo->lL)][mu][col].re = 1.;


	for (i = 0; i < nsmear; i++)
	{
		if (qcd_momgaussIteration3d(vec, uAPE, alpha, x_src[0], mom, keci)) //注意动量的数据类型
		{
			fprintf(stderr, "process %i: Error while smearing!\n", geo->myid);
			exit(EXIT_FAILURE);
		}
	}

	//qcd_writeVectorLime(vec_names[is], qcd_SOURCE_LIME, &vec);
	if (geo->myid == 0) printf(" Done vector: %4d / %4d\n", is, numb_sources);

}

void WallSource3D(qcd_vector* vec, qcd_int_4 mom[3], qcd_int_4 t,qcd_int_4 is)
{
	qcd_int_4 lt,lx, ly, lz;
	qcd_int_4 x, y, z;
	qcd_int_8 v;
	qcd_int_4 mu = is / 3;
	qcd_int_4 col = is % 3;
	qcd_real_8 tmp;
	qcd_zeroVector(vec);
	
	lt = t - vec->geo->Pos[0] * vec->geo->lL[0];

	if (lt >= 0 && lt < vec->geo->lL[0])
	{
		for (lx = 0; lx < vec->geo->lL[1]; lx++)
		for (ly = 0; ly < vec->geo->lL[2]; ly++)
		for (lz = 0; lz < vec->geo->lL[3]; lz++)
		{
			v = qcd_LEXIC(lt, lx, ly, lz, vec->geo->lL);
			x = lx + vec->geo->Pos[1] * vec->geo->lL[1];
			y = ly + vec->geo->Pos[2] * vec->geo->lL[2];
			z = lz + vec->geo->Pos[3] * vec->geo->lL[3];
			tmp = (((double)mom[0] * x) / vec->geo->L[1] + ((double)mom[1] * y) / vec->geo->L[2] + ((double)mom[2] * z) / vec->geo->L[3]) * 2 * M_PI;
			vec->D[v][mu][col] = (qcd_complex_16){ cos(tmp), sin(tmp) };
		}
	}

}


void WallSource4D(qcd_vector* vec, qcd_int_4 mom[4], qcd_int_4 is)//重整化off shell 时需要,四维的傅里叶变换
{
	qcd_int_4 lt, lx, ly, lz;
	qcd_int_4 t,x, y, z;
	qcd_int_8 v;
	qcd_int_4 mu = is / 3;
	qcd_int_4 col = is % 3;
	qcd_real_8 tmp;
	qcd_complex_16 phase_factor1, phase_factor2;
	qcd_zeroVector(vec);

	for (lt = 0; lt < vec->geo->lL[0]; lt++)
	for (lx = 0; lx < vec->geo->lL[1]; lx++)
	for (ly = 0; ly < vec->geo->lL[2]; ly++)
	for (lz = 0; lz < vec->geo->lL[3]; lz++)
	{
		v = qcd_LEXIC(lt, lx, ly, lz, vec->geo->lL);
		t = lt + vec->geo->Pos[0] * vec->geo->lL[0];
		x = lx + vec->geo->Pos[1] * vec->geo->lL[1];
		y = ly + vec->geo->Pos[2] * vec->geo->lL[2];
		z = lz + vec->geo->Pos[3] * vec->geo->lL[3];
		tmp = ( ((double) mom[0]*t) /vec->geo->L[0]+ ((double)mom[1] * x) / vec->geo->L[1] + \
			((double)mom[2] * y) / vec->geo->L[2] + ((double)mom[3] * z) / vec->geo->L[3]) * 2 * M_PI;
		
		phase_factor1 = (qcd_complex_16){ cos(tmp), sin(tmp) };
		//phase_factor2 = (qcd_complex_16){ cos(M_PI * t / vec->geo->L[0]), -sin(M_PI * t / vec->geo->L[0]) };
		phase_factor2 = (qcd_complex_16){ 1,0 };
		vec->D[v][mu][col] =qcd_CMUL(phase_factor1,phase_factor2) ;
	}
	

}




void mksource_stochasticZ4(qcd_vector *vec,qcd_int_4 t,qcd_int_4 is)
{
	int mu, col;

	qcd_int_4 lt,lx, ly, lz, x, y, z;
	qcd_real_8 tmp;
	qcd_int_8 v;
	qcd_complex_16 momfactor;
	qcd_complex_16 randompoint;

	mu = is/3;
	col = is%3;

	qcd_zeroVector(vec);

	lt = t - vec->geo->Pos[0] * vec->geo->lL[0];

	if (lt >= 0 && lt < vec->geo->lL[0])
	{
		for (lx = 0; lx < vec->geo->lL[1]; lx++)
		for (ly = 0; ly < vec->geo->lL[2]; ly++)
		for (lz = 0; lz < vec->geo->lL[3]; lz++)
		{
			v = qcd_LEXIC(lt, lx, ly, lz, vec->geo->lL);
			randompoint = (qcd_complex_16) { pow(-1, rand()), pow(-1, rand()) };
			vec->D[v][mu][col] = qcd_CSCALE( randompoint,sqrt(2)/2);
		}
	}

}


void stochasticZ4(qcd_vector* des, qcd_vector* vec, qcd_int_4 mom[3], qcd_int_4 t,qcd_int_4 dis,qcd_int_4 vis)
{
	int mu, col;

	qcd_int_4 lt, lx, ly, lz, x, y, z;
	qcd_real_8 tmp;
	qcd_int_8 v;
	qcd_complex_16 momfactor;
	qcd_complex_16 randompoint;

	mu = dis /3;
	col = dis%3;

	qcd_zeroVector(des);

	lt = t - vec->geo->Pos[0] * vec->geo->lL[0];

	if (lt >= 0 && lt < vec->geo->lL[0])
	{
		for (lx = 0; lx < vec->geo->lL[1]; lx++)
		for (ly = 0; ly < vec->geo->lL[2]; ly++)
		for (lz = 0; lz < vec->geo->lL[3]; lz++)
		{
			v = qcd_LEXIC(lt, lx, ly, lz, vec->geo->lL);
			x = lx + vec->geo->Pos[1] * vec->geo->lL[1];
			y = ly + vec->geo->Pos[2] * vec->geo->lL[2];
			z = lz + vec->geo->Pos[3] * vec->geo->lL[3];
			tmp = (((double)mom[0] * x) / vec->geo->L[1] + ((double)mom[1] * y) / vec->geo->L[2] + ((double)mom[2] * z) / vec->geo->L[3]) * 2 * M_PI;
			momfactor = (qcd_complex_16){ cos(tmp), sin(tmp) };
			
			des->D[v][mu][col] = qcd_CMUL(momfactor, vec->D[v][vis/3][vis%3]); //vis配合mksource_stochasticZ4 使用
		}
	}

}