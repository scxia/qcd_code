#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <qcd.h>
#include <qcd_projectors.h>

void twop_proton_field_selection(qcd_propagator * uprop, qcd_propagator* dprop, qcd_geometry *geo, qcd_int_4 x_src[4], qcd_int_2 t_start, \
	qcd_int_2 t_stop, qcd_int_4 mom[3], FILE *fp_corr_p,qcd_int_4 (*random_list)[4],qcd_int_4 random_num)
{
	qcd_int_4 ctr, ctr2, cc1, cc2, lx, ly, lz, x, y, z;
	qcd_int_4 c1, c2, c3, c1p, c2p, c3p;
	qcd_int_4 mu, nu, ku, lu, gi, al, be;
	qcd_int_4 random_index;
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
	for (t = t_start; t <= t_stop; t=t+2)
	{
		lt = ((t + x_src[0]) % geo->L[0]) - geo->Pos[0] * geo->lL[0];

		//if (geo->myid == 0) printf("t=%i\n", t);
		for (v3 = 0; v3 < geo->lV3; v3++)   //set blocks to zero
			block[v3] = (qcd_complex_16) { 0, 0 };

		if (lt >= 0 && lt < geo->lL[0])  //inside the local lattice, otherwise nothing to calculate
		{


			for (lx = 0; lx < geo->lL[1]; lx++)
				for (ly = 0; ly < geo->lL[2]; ly++)
				for (lz = 0; lz < geo->lL[3]; lz++)
				{



			      v3 = qcd_LEXIC0(lx, ly, lz, geo->lL);
			      v = qcd_LEXIC(lt, lx, ly, lz, geo->lL);


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



									block[v3] = qcd_CSUB(block[v3]
										, qcd_CMUL(qcd_ONE_PLUS_GAMMA[0][al][be]
											, qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(cgicgi_val[ctr2]
												, uprop->D[v][mu][ku][c1][c1p])
												, dprop->D[v][nu][lu][c2][c2p])
												, uprop->D[v][be][al][c3][c3p])
												, qcd_SGN_EPS[cc1] * qcd_SGN_EPS[cc2])));

									block[v3] = qcd_CSUB(block[v3]
										, qcd_CMUL(qcd_ONE_PLUS_GAMMA[0][al][be]
											, qcd_CSCALE(qcd_CMUL(qcd_CMUL(qcd_CMUL(cgicgi_val[ctr2]
												, uprop->D[v][mu][al][c1][c1p])
												, dprop->D[v][nu][lu][c2][c2p])
												, uprop->D[v][be][ku][c3][c3p])
												, qcd_SGN_EPS[cc1] * qcd_SGN_EPS[cc2])));


								}
							}			//space loop
						}//color2 loop    
					}			//color1 loop
				}//nonvanishing cg5cg5 loop
			//nonvanishing projector condition
		}
		 //Fourier transform time-slice

		//if (geo->myid == 0) printf("now fourier trans\n");
			//for (j = 0; j < momi; j++)

		for (random_index = 0; random_index < random_num; random_index++)
		{
			x = random_list[((t + x_src[0]) % geo->L[0])/2 * random_num + random_index][1];
			y = random_list[((t + x_src[0]) % geo->L[0])/2 * random_num + random_index][2];
			z = random_list[((t + x_src[0]) % geo->L[0])/2 * random_num + random_index][3];
			lx = x - geo->Pos[1] * geo->lL[1];
			ly = y - geo->Pos[2] * geo->lL[2];
			lz = z - geo->Pos[3] * geo->lL[3];
			
			//{
			if (geo->myid == 0)
			{
				fprintf(fp_corr_p, "%i %+i %+i %+i ", t, mom[0], mom[1], mom[2]);
			}
			corr = (qcd_complex_16) { 0, 0 };

			if ((lx >= 0) && (lx < geo->lL[1]) && (ly >= 0) && (ly < geo->lL[2]) && (lz >= 0) && (lz < geo->lL[3]))
			{
				v3 = qcd_LEXIC0(lx, ly, lz, geo->lL);
				x = lx + geo->Pos[1] * geo->lL[1] - x_src[1];
				y = ly + geo->Pos[2] * geo->lL[2] - x_src[2];
				z = lz + geo->Pos[3] * geo->lL[3] - x_src[3];
				tmp = (((double)mom[0] * x) / geo->L[1] + ((double)mom[1] * y) / geo->L[2] + ((double)mom[2] * z) / geo->L[3]) * 2 * M_PI;
				C2 = (qcd_complex_16) { cos(tmp), -sin(tmp) }; //TABULATE FOR LARGE SPEEDUP!!!
				corr = qcd_CADD(corr, qcd_CMUL(block[v3], C2));

			}
			//if (geo->myid == 0) printf("now fourier trans finish \n");
			//printf("process %i: corr = %f %+fi\n",myid,0.5*corr.re,0.5*corr.im);
			MPI_Reduce(&(corr.re), &(corr2.re), 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			if (geo->myid == 0)
			{
				fprintf(fp_corr_p, "%+e %+e\n", corr2.re*0.25, corr2.im*0.25);  //投影算符中漏的1/4
			}
		}
		//if (geo->myid == 0) printf("one corr finish\n");
	}

	//end lt inside local block condition

	free(block);


}






void mksource_su_proton_field_selection(qcd_propagator *suprop, qcd_propagator *uprop, qcd_propagator *dprop, qcd_propagator *upropsink, qcd_propagator *dpropsink, qcd_gaugeField *uAPE, qcd_geometry *geo, \
	qcd_int_4 x_src[4], qcd_real_8 alpha, qcd_uint_4 nsmear, qcd_real_8 keci, qcd_int_4 mom[3], qcd_int_4 tsink,qcd_int_4 randompoint[4])
{
	qcd_int_4 mu, nu, ku, lu, be, al, ta;
	qcd_int_4 c1, c2, c3, c1p, c2p, c3p, cc1, cc2, a, b;
	qcd_int_4 ctr, ctr2, row, col, i, gi, v;
	qcd_int_4 x, y, z, lx, ly, lz, t, lt;
	//qcd_int_4 numb_sources = 12;

	qcd_real_8 tmp;
	qcd_complex_16 C, C2;
	qcd_complex_16 imag = { 0, 1 };

	qcd_int_2 cgicgi_ind[16 * 16 * 16][6];
	qcd_complex_16 cgicgi_val[16 * 16 * 16];
	qcd_complex_16 sou[4][4][3][3] ;
	qcd_complex_16 sou2[4][4][3][3];

	memset(sou, 0, 4 * 4 * 3 * 3 * sizeof(qcd_complex_16));
	memset(sou2, 0, 4 * 4 * 3 * 3 * sizeof(qcd_complex_16));

	//qcd_propagator soupro;

	//numb_sources = 12;


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
	//qcd_initPropagator(&soupro, geo);
	//qcd_zeroPropagator(&soupro);
	//qcd_zeroVector(souvec);
	t = tsink;
	lt = ((t + x_src[0]) % geo->L[0]) - geo->Pos[0] * geo->lL[0];
	lx = randompoint[1] - geo->Pos[1] * geo->lL[1];
	ly = randompoint[2] - geo->Pos[2] * geo->lL[2];
	lz = randompoint[3] - geo->Pos[3] * geo->lL[3];
	
	if (lt >= 0 && lt < geo->lL[0] &&lx>=0 && lx<geo->lL[1]&& ly >= 0 && ly < geo->lL[2] && lz >= 0 && lz < geo->lL[3])  //inside the local lattice, otherwise nothing to calculate
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
					//for (ta = 0; ta < 4; ta++)
					//{
						//if (qcd_NORM(qcd_CADD(qcd_GAMMA[5][ta][mu], qcd_CMUL(imag, PROJECTOR[0][ta][mu]))) > 1e-3)
						//{
							//for (lx = 0; lx < geo->lL[1]; lx++)
								//for (ly = 0; ly < geo->lL[2]; ly++)
									//for (lz = 0; lz < geo->lL[3]; lz++)
									//{
										//v3 = qcd_LEXIC0(lx, ly, lz, geo.lL);
										v = qcd_LEXIC(lt, lx, ly, lz, geo->lL);

										sou[mu][lu][c1][c3p] = qcd_CADD(sou[mu][lu][c1][c3p], 
												qcd_CSCALE(qcd_CMUL(qcd_CMUL(cgicgi_val[ctr2]
													, dprop->D[v][nu][ku][c2][c2p])
													, uprop->D[v][al][be][c3][c1p])
													, qcd_SGN_EPS[cc1] * qcd_SGN_EPS[cc2]));

									
									//}
					//	}
					//}//space loop
				}//color2 loop    
			}//color1 loop
		}//nonvanishing cg5cg5 loop
	//nonvanishing projector condition

	//Fourier transform time-slice


	//	printf("sou re +%e im +%e \n", sou[0][0][0][0].re, sou[0][0][0][0].im);
	//	printf("sou re +%e im +%e \n", sou[2][2][2][2].re, sou[2][2][2][2].im);

	}
	   
	MPI_Allreduce(&(sou[0][0][0][0].re), &(sou2[0][0][0][0].re),4*4*3*3*2,MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);

	for (v = 0; v < geo->lV; v++)
	for (mu = 0; mu < 4; mu++)
	for (nu = 0; nu < 4; nu++)
	for (al = 0; al < 4; al++)
	for (c1 = 0; c1 < 3; c1++)
	for (c1p = 0; c1p < 3; c1p++)
	for (c2 = 0; c2 < 3; c2++)
	{
		suprop->D[v][mu][al][c1][c2] = qcd_CADD(suprop->D[v][mu][al][c1][c2],
			qcd_CMUL(qcd_CMUL(qcd_CONJ(dpropsink->D[v][mu][nu][c1][c1p]),
				qcd_CMUL(qcd_GAMMA[5][mu][mu], qcd_GAMMA[5][nu][nu])), sou2[nu][al][c1p][c2]));
	}


		

}







void mksource_sd_proton_field_selection(qcd_propagator *sdprop, qcd_propagator *uprop, qcd_propagator *dprop, qcd_propagator *upropsink, qcd_propagator *dpropsink, qcd_gaugeField *uAPE, qcd_geometry *geo, \
	qcd_int_4 x_src[4], qcd_real_8 alpha, qcd_uint_4 nsmear, qcd_real_8 keci, qcd_int_4 mom[3], qcd_int_4 tsink, qcd_int_4 randompoint[4])
{

	qcd_int_4 mu, nu, ku, lu, be, al, ta;
	qcd_int_4 c1, c2, c3, c1p, c2p, c3p, cc1, cc2, a, b;
	qcd_int_4 ctr, ctr2, row, col, i, gi, v;
	qcd_int_4 x, y, z, lx, ly, lz, t, lt;
	//qcd_int_4 numb_sources = 12;

	qcd_real_8 tmp;
	qcd_complex_16 C, C2;
	qcd_complex_16 imag = { 0, 1 };

	qcd_int_2 cgicgi_ind[16 * 16 * 16][6];
	qcd_complex_16 cgicgi_val[16 * 16 * 16];


	qcd_complex_16 sou[4][4][3][3];
	qcd_complex_16 sou2[4][4][3][3];

	memset(sou, 0, 4 * 4 * 3 * 3*sizeof(qcd_complex_16) );
	memset(sou2, 0, 4 * 4 * 3 * 3 * sizeof(qcd_complex_16));

	//numb_sources = 12;

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
	//qcd_initPropagator(&soupro, geo);
	//qcd_zeroPropagator(&soupro);
	//qcd_zeroVector(souvec);
	t = tsink;
	lt = ((t + x_src[0]) % geo->L[0]) - geo->Pos[0] * geo->lL[0];
	lx = randompoint[1] - geo->Pos[1] * geo->lL[1];
	ly = randompoint[2] - geo->Pos[2] * geo->lL[2];
	lz = randompoint[3] - geo->Pos[3] * geo->lL[3];

	if (lt >= 0 && lt < geo->lL[0]&& lx >= 0 && lx < geo->lL[1] && ly >= 0 && ly < geo->lL[2] && lz >= 0 && lz < geo->lL[3])  //inside the local lattice, otherwise nothing to calculate
	{
		//if (geo->myid == 0) printf("t=%i\n", t);

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
					/*for (ta = 0; ta < 4; ta++)
					{
						if (qcd_NORM(qcd_CSUB(qcd_GAMMA[5][ta][nu], qcd_CMUL(imag, PROJECTOR[0][ta][nu]))) > 1e-3)
						{
							for (lx = 0; lx < geo->lL[1]; lx++)
								for (ly = 0; ly < geo->lL[2]; ly++)
									for (lz = 0; lz < geo->lL[3]; lz++)
									{*/
										//v3 = qcd_LEXIC0(lx, ly, lz, geo.lL);
										v = qcd_LEXIC(lt, lx, ly, lz, geo->lL);

										sou[nu][ku][c2][c2p] = qcd_CADD(sou[nu][ku][c2][c2p],
											 qcd_CSCALE(qcd_CMUL(cgicgi_val[ctr2],
												qcd_CSUB(qcd_CMUL(uprop->D[v][mu][lu][c1][c3p], uprop->D[v][al][be][c3][c1p]), qcd_CMUL(uprop->D[v][mu][be][c1][c1p], uprop->D[v][al][lu][c3][c3p])))
												, qcd_SGN_EPS[cc1] * qcd_SGN_EPS[cc2]));


									

									

								//	}
					//	}
					//}//space loop
				}//color2 loop    
			}//color1 loop
		}//nonvanishing cg5cg5 loop
	//nonvanishing projector condition

		//printf("sou re +%e im +%e \n", sou[0][0][0][0].re, sou[0][0][0][0].im);
		//printf("sou re +%e im +%e \n", sou[2][2][2][2].re, sou[2][2][2][2].im);
		

	}

	MPI_Allreduce(&(sou[0][0][0][0].re), &(sou2[0][0][0][0].re), 4 * 4 * 3 * 3 * 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	for (v = 0; v < geo->lV; v++)
		for (mu = 0; mu < 4; mu++)
			for (nu = 0; nu < 4; nu++)
				for (al = 0; al < 4; al++)
					for (c1 = 0; c1 < 3; c1++)
						for (c1p = 0; c1p < 3; c1p++)
							for (c2 = 0; c2 < 3; c2++)
							{
								sdprop->D[v][mu][al][c1][c2] = qcd_CADD(sdprop->D[v][mu][al][c1][c2],
									qcd_CMUL(qcd_CMUL(qcd_CONJ(upropsink->D[v][mu][nu][c1][c1p]),
										qcd_CMUL(qcd_GAMMA[5][mu][mu], qcd_GAMMA[5][nu][nu])), sou2[nu][al][c1p][c2]));
							}





	
}




void axial_charge_field_selection(qcd_propagator *uprop, qcd_propagator *dprop, qcd_propagator*usprop, qcd_propagator *dsprop, \
	qcd_geometry *geo, qcd_int_4 x_src[4], qcd_int_4 t_start, qcd_int_4 t_stop, FILE *fp_corr_p, qcd_int_4 mom[3],qcd_int_4 (*random_list)[4],qcd_int_4 random_num,qcd_int_4 nsou,qcd_int_4 nsink)
{
	qcd_int_4 mu, ku, nu, lu;
	qcd_int_4 x, y, z, t, lx, ly, lz, lt;
	qcd_int_4 cc1, cc2, cc3, wl;
	qcd_int_8 v, v3, vm, vp;
	qcd_int_4 random_index;
	qcd_int_8 ctrm, ctrm2;

	qcd_real_8 tmp;
	qcd_complex_16 C, C2;
	qcd_complex_16 corr, corr2;
	qcd_int_2 cg5cg5m_ind[16 * 16][4];
	qcd_complex_16 cg5cg5m_val[16 * 16];
	qcd_complex_16 imag = { 0,1 };


	qcd_complex_16 *block;


	block = (qcd_complex_16*)malloc(geo->lV3 * sizeof(qcd_complex_16));


	ctrm = 0;
	for (mu = 0; mu < 4; mu++)
			for (nu = 0; nu < 4; nu++)
			{
				C = qcd_CMUL(imag,qcd_GAMMAG5[3][mu][nu]);   
				if (qcd_NORM(C) > 1e-3)
				{
					cg5cg5m_val[ctrm].re = C.re;
					cg5cg5m_val[ctrm].im = C.im;
					cg5cg5m_ind[ctrm][0] = mu;
					cg5cg5m_ind[ctrm][1] = nu;
					//cg5cg5m_ind[ctrm][2] = ku;
					ctrm++;
				}
			}

	int selection_num[15] = { 13824,1728,864,432,216,108,64,54,32,27,16,8,4,2,1 };
	int ss = 0;

		for (t = t_start; t <= t_stop; t++)
		{
			lt = ((t + x_src[0]) % geo->L[0]) - geo->Pos[0] * geo->lL[0];

			for (v3 = 0; v3 < geo->lV3; v3++)   //set blocks to zero
				block[v3] = (qcd_complex_16) { 0, 0 };

			if (lt >= 0 && lt < geo->lL[0])  //inside the local lattice, otherwise nothing to calculate
			{

				//if (geo->myid == 0) printf("time %i\n", t);
				
				
					for (lx = 0; lx < geo->lL[1]; lx++)
						for (ly = 0; ly < geo->lL[2]; ly++)
							for (lz = 0; lz < geo->lL[3]; lz++)
							{

								v = qcd_LEXIC(lt, lx, ly, lz, geo->lL);
								v3 = qcd_LEXIC0(lx, ly, lz, geo->lL);



								for (lu = 0; lu < 4; lu++)
								{
									for (cc1 = 0; cc1 < 3; cc1++)
									{
										for (cc2 = 0; cc2 < 3; cc2++)
										{

											for (ctrm2 = 0; ctrm2 < ctrm; ctrm2++)
											{
												mu = cg5cg5m_ind[ctrm2][0];
												nu = cg5cg5m_ind[ctrm2][1];
												//ku = cg5cg5m_ind[ctrm2][2];

												block[v3] = qcd_CADD(block[v3]
													, qcd_CMUL(qcd_CMUL(cg5cg5m_val[ctrm2]
														, uprop->D[v][nu][lu][cc1][cc2])
														, usprop->D[v][mu][lu][cc1][cc2]));

												block[v3] = qcd_CSUB(block[v3]
													, qcd_CMUL(qcd_CMUL(cg5cg5m_val[ctrm2]
														, dprop->D[v][nu][lu][cc1][cc2])
														, dsprop->D[v][mu][lu][cc1][cc2]));

											}

										}//space loop
									}//color2 loop    
								}//color1 loop
								/*x = lx + geo->Pos[1] * geo->lL[1];
								y = ly + geo->Pos[2] * geo->lL[2];
								z = lz + geo->Pos[3] * geo->lL[3];
								if (t==6&&x == 12 && y == 13 && z == 13)
								{
									printf("blockv3 %d %d %d %e %e\n", x, y, z, block[v3].re,block[v3].im);
								}*/


							}


			}
			//Fourier transform time-slice

			for (ss = 0; ss < 15; ss++)
			{

				if (geo->myid == 0)
				{
					// fprintf(fp_corr_p, "%+i %+i %+i %+i %+i", wl, mom[j][0], mom[j][1], mom[j][2],wl);
					fprintf(fp_corr_p, "%2d %2d %2d %5d ", nsou, nsink, t, selection_num[ss]);

				}
				corr = (qcd_complex_16) { 0, 0 };


				if (ss == 0)
				{
					for (lx = 0; lx < geo->lL[1]; lx++)
						for (ly = 0; ly < geo->lL[2]; ly++)
							for (lz = 0; lz < geo->lL[3]; lz++)
							{
								/*x = lx + geo->Pos[1] * geo->lL[1];
								y = ly + geo->Pos[2] * geo->lL[2];
								z = lz + geo->Pos[3] * geo->lL[3];
								if (t == 6 && x == 12 && y == 13 && z == 13)
								{
									printf("Fourier blockv3 %d %d %d %e %e\n", x, y, z, block[v3].re, block[v3].im);
								}*/
								v3 = qcd_LEXIC0(lx, ly, lz, geo->lL);
								x = lx + geo->Pos[1] * geo->lL[1] - x_src[1];
								y = ly + geo->Pos[2] * geo->lL[2] - x_src[2];
								z = lz + geo->Pos[3] * geo->lL[3] - x_src[3];
								tmp = (((double)mom[0] * x) / geo->L[1] + ((double)mom[1] * y) / geo->L[2] + ((double)mom[2] * z) / geo->L[3]) * 2 * M_PI;
								C2 = (qcd_complex_16) { cos(tmp), -sin(tmp) }; //TABULATE FOR LARGE SPEEDUP!!!
								corr = qcd_CADD(corr, qcd_CMUL(block[v3], C2));
							}
				}
				else
				{
					for (random_index = 0; random_index < selection_num[ss]; random_index++)
					{
						x = random_list[((t + x_src[0]) % geo->L[0]) * random_num + random_index][1];
						y = random_list[((t + x_src[0]) % geo->L[0])* random_num + random_index][2];
						z = random_list[((t + x_src[0]) % geo->L[0]) * random_num + random_index][3];
						lx = x - geo->Pos[1] * geo->lL[1];
						ly = y - geo->Pos[2] * geo->lL[2];
						lz = z - geo->Pos[3] * geo->lL[3];
						if ((lx >= 0) && (lx < geo->lL[1]) && (ly >= 0) && (ly < geo->lL[2]) && (lz >= 0) && (lz < geo->lL[3]))
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
				}

				//MPI_Allreduce(&(corr.re), &(corr2.re), 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);


				MPI_Reduce(&(corr.re), &(corr2.re), 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

				if (geo->myid == 0)
				{
					fprintf(fp_corr_p, "%+e %+e\n", corr2.re, corr2.im);
				}
			}



		}
	

	free(block);

}