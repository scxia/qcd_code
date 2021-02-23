#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <qcd.h>
#include <qcd_projectors.h>
#include<time.h>



void axial_charge(qcd_propagator *uprop, qcd_propagator *dprop, qcd_propagator*usprop, qcd_propagator *dsprop, \
	 qcd_geometry *geo, qcd_int_4 x_src[4], qcd_int_4 t_start, qcd_int_4 t_stop, \
	 FILE *fp_corr_p, qcd_int_4 mom[3])
{
	qcd_int_4 mu, ku, nu, lu;
	qcd_int_4 x, y, z, t, lx, ly, lz, lt;
	qcd_int_4 cc1, cc2, cc3, wl;
	qcd_int_8 v, v3, vm, vp;
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
		for (ku = 0; ku < 4; ku++)
			for (nu = 0; nu < 4; nu++)
			{
				C = qcd_CMUL(qcd_GAMMAG5[3][mu][nu], qcd_GAMMA[5][mu][ku]);   //g5 from the sequential source -> D_u =g5 * D_d^daggar * g5
				if (qcd_NORM(C) > 1e-3)
				{
					cg5cg5m_val[ctrm].re = C.re;
					cg5cg5m_val[ctrm].im = C.im;
					cg5cg5m_ind[ctrm][0] = mu;
					cg5cg5m_ind[ctrm][1] = nu;
					cg5cg5m_ind[ctrm][2] = ku;
					ctrm++;
				}
			}



	for (t = t_start; t <= t_stop; t++)
	{
		lt = ((t + x_src[0]) % geo->L[0]) - geo->Pos[0] * geo->lL[0];

		for (v3 = 0; v3 < geo->lV3; v3++)   //set blocks to zero
			block[v3] = (qcd_complex_16) { 0, 0 };

		if (lt >= 0 && lt < geo->lL[0])  //inside the local lattice, otherwise nothing to calculate
		{

			if (geo->myid == 0) printf("time %i\n", t);

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
										ku = cg5cg5m_ind[ctrm2][2];




										block[v3] = qcd_CADD(block[v3]
											, qcd_CMUL(qcd_CMUL(cg5cg5m_val[ctrm2]
												, uprop->D[v][nu][lu][cc1][cc2])
												, qcd_CONJ(qcd_CMUL(qcd_CSUB(PROJECTOR[0][ku][ku]
													, qcd_CMUL(imag, qcd_GAMMA[5][ku][ku])), usprop->D[v][ku][lu][cc1][cc2]))));

										block[v3] = qcd_CSUB(block[v3]
											, qcd_CMUL(qcd_CMUL(cg5cg5m_val[ctrm2]
												, dprop->D[v][nu][lu][cc1][cc2])
												, qcd_CONJ(qcd_CMUL(qcd_CADD(PROJECTOR[0][ku][ku]
													, qcd_CMUL(imag, qcd_GAMMA[5][ku][ku])), dsprop->D[v][ku][lu][cc1][cc2]))));

									}

								}//space loop
							}//color2 loop    
						}//color1 loop


					}
			//nonvanishing cg5cg5 loop
			//nonvanishing projector condition

		}
		//Fourier transform time-slice

		
		
			if (geo->myid == 0)
			{
				// fprintf(fp_corr_p, "%+i %+i %+i %+i %+i", wl, mom[j][0], mom[j][1], mom[j][2],wl);
				fprintf(fp_corr_p, "%d ", t);

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
			MPI_Reduce(&(corr.re), &(corr2.re), 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

			if (geo->myid == 0)
			{
				fprintf(fp_corr_p, "%+e %+e\n", corr2.re*0.5, corr2.im*0.5);
			}
		


	}

	free(block);

}






void threepx(qcd_propagator *uprop, qcd_propagator *dprop, qcd_propagator*usprop, qcd_propagator *dsprop, \
	qcd_gaugeField *ustout, qcd_geometry *geo, qcd_int_4 x_src[4], qcd_int_4 t_start, qcd_int_4 t_stop, \
	qcd_int_4 w_length, FILE *fp_corr_p, qcd_int_4 mom[3])
{

	qcd_int_4 mu, ku, nu, lu;
	qcd_int_4 x, y, z, t, lx, ly, lz, lt;
	qcd_int_4 cc1, cc2, cc3, wl;
	qcd_int_8 v, v3, vm, vp;
	qcd_int_8 ctrm, ctrm2;

	qcd_real_8 tmp;
	qcd_complex_16 C, C2;
	qcd_complex_16 corr, corr2;
	qcd_int_2 cg5cg5m_ind[16 * 16][4];
	qcd_complex_16 cg5cg5m_val[16 * 16];
	qcd_complex_16 imag = { 0,1 };


	qcd_complex_16 *block;
	qcd_complex_16(*wilsonline)[3][3];
	qcd_complex_16 tempp[3][3];


	block = (qcd_complex_16*)malloc(geo->lV3 * sizeof(qcd_complex_16));
	wilsonline = malloc(geo->lV * 9 * sizeof(qcd_complex_16));

	if (wilsonline == NULL)
	{
		fprintf(stderr, "process %i: Error in wilsonline, Out of memory\n", geo->myid); exit(EXIT_FAILURE);
	}
	else
	{
		if (geo->myid == 0) printf("the wilsonline init sucessfully\n");
	}


	ctrm = 0;
	for (mu = 0; mu < 4; mu++)
		for (ku = 0; ku < 4; ku++)
			for (nu = 0; nu < 4; nu++)
			{
				C = qcd_CMUL(qcd_GAMMA[0][mu][nu], qcd_GAMMA[5][mu][ku]);
				if (qcd_NORM(C) > 1e-3)
				{
					cg5cg5m_val[ctrm].re = C.re;
					cg5cg5m_val[ctrm].im = C.im;
					cg5cg5m_ind[ctrm][0] = mu;
					cg5cg5m_ind[ctrm][1] = nu;
					cg5cg5m_ind[ctrm][2] = ku;
					ctrm++;
				}
			}

	for (v3 = 0; v3 < geo->lV; v3++)
	{
		qcd_unit3x3(wilsonline[v3]);
	}

	for (wl = 0; wl <= w_length; wl++)
	{
		for (t = t_start; t <= t_stop; t++)
		{
			lt = ((t + x_src[0]) % geo->L[0]) - geo->Pos[0] * geo->lL[0];

			for (v3 = 0; v3 < geo->lV3; v3++)   //set blocks to zero
				block[v3] = (qcd_complex_16) { 0, 0 };

			if (geo->myid == 0) printf(" calculating wilesonline %i time %i\n", wl, t);

			if (lt >= 0 && lt < geo->lL[0])  //inside the local lattice, otherwise nothing to calculate
			{

				for (lx = 0; lx < geo->lL[1]; lx++)
					for (ly = 0; ly < geo->lL[2]; ly++)
						for (lz = 0; lz < geo->lL[3]; lz++)
						{
							vp = qcd_LEXIC(lt, (lx + wl - 1) % geo->lL[1], ly, lz, geo->lL);
							vm = qcd_LEXIC(lt, (lx + wl) % geo->lL[1], ly, lz, geo->lL);
							v = qcd_LEXIC(lt, lx, ly, lz, geo->lL);
							v3 = qcd_LEXIC0(lx, ly, lz, geo->lL);
							if (wl != 0)
							{
								qcd_mul3x3(tempp, wilsonline[v], ustout->D[vp][1]);
								qcd_copy3x3(wilsonline[v], tempp);
							}
							for (lu = 0; lu < 4; lu++)
							{
								for (cc1 = 0; cc1 < 3; cc1++)
								{
									for (cc2 = 0; cc2 < 3; cc2++)
									{
										for (cc3 = 0; cc3 < 3; cc3++)
										{
											for (ctrm2 = 0; ctrm2 < ctrm; ctrm2++)
											{
												mu = cg5cg5m_ind[ctrm2][0];
												nu = cg5cg5m_ind[ctrm2][1];
												ku = cg5cg5m_ind[ctrm2][2];




												block[v3] = qcd_CADD(block[v3]
													, qcd_CMUL(wilsonline[v][cc1][cc3], qcd_CMUL(qcd_CMUL(cg5cg5m_val[ctrm2]
														, uprop->D[vm][nu][lu][cc3][cc2])
														, qcd_CONJ(qcd_CMUL(qcd_CSUB(PROJECTOR[0][ku][ku]
															, qcd_CMUL(imag, qcd_GAMMA[5][ku][ku])), usprop->D[v][ku][lu][cc1][cc2])))));

												block[v3] = qcd_CSUB(block[v3]
													, qcd_CMUL(wilsonline[v][cc1][cc3], qcd_CMUL(qcd_CMUL(cg5cg5m_val[ctrm2]
														, dprop->D[vm][nu][lu][cc3][cc2])
														, qcd_CONJ(qcd_CMUL(qcd_CADD(PROJECTOR[0][ku][ku]
															, qcd_CMUL(imag, qcd_GAMMA[5][ku][ku])), dsprop->D[v][ku][lu][cc1][cc2])))));

											}
										}
									}//space loop
								}//color2 loop    
							}//color1 loop

						}
			}



			if (geo->myid == 0)
			{
				// fprintf(fp_corr_p, "%+i %+i %+i %+i %+i", wl, mom[j][0], mom[j][1], mom[j][2],wl);
				fprintf(fp_corr_p, "%d  %d ", wl, t);

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

			MPI_Reduce(&(corr.re), &(corr2.re), 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			if (geo->myid == 0)
			{
				fprintf(fp_corr_p, "%+e %+e\n", corr2.re*0.5, corr2.im*0.5);
			}



		}//end lt inside local block condition

	}//end t-loop   

	for (v3 = 0; v3 < geo->lV; v3++)
	{
		qcd_unit3x3(wilsonline[v3]);
	}

	for (wl = -1; wl >= -w_length; wl--)
	{
		for (t = t_start; t <= t_stop; t++)
		{
			lt = ((t + x_src[0]) % geo->L[0]) - geo->Pos[0] * geo->lL[0];

			for (v3 = 0; v3 < geo->lV3; v3++)   //set blocks to zero
				block[v3] = (qcd_complex_16) { 0, 0 };

			if (geo->myid == 0) printf("calculating wilesonline %i time %i\n", wl, t);

			if (lt >= 0 && lt < geo->lL[0])  //inside the local lattice, otherwise nothing to calculate
			{

				for (lx = 0; lx < geo->lL[1]; lx++)
					for (ly = 0; ly < geo->lL[2]; ly++)
						for (lz = 0; lz < geo->lL[3]; lz++)
						{
							vp = qcd_LEXIC(lt, (lx + wl + geo->lL[1]) % geo->lL[1], ly, lz, geo->lL);
							v = qcd_LEXIC(lt, lx, ly, lz, geo->lL);
							v3 = qcd_LEXIC0(lx, ly, lz, geo->lL);
							qcd_mulAdjoint3x3(tempp, wilsonline[v], ustout->D[vp][1]);
							qcd_copy3x3(wilsonline[v], tempp);

							for (lu = 0; lu < 4; lu++)
							{
								for (cc1 = 0; cc1 < 3; cc1++)
								{
									for (cc2 = 0; cc2 < 3; cc2++)
									{
										for (cc3 = 0; cc3 < 3; cc3++)
										{
											for (ctrm2 = 0; ctrm2 < ctrm; ctrm2++)
											{
												mu = cg5cg5m_ind[ctrm2][0];
												nu = cg5cg5m_ind[ctrm2][1];
												ku = cg5cg5m_ind[ctrm2][2];


												block[v3] = qcd_CADD(block[v3]
													, qcd_CMUL(wilsonline[v][cc1][cc3], qcd_CMUL(qcd_CMUL(cg5cg5m_val[ctrm2]
														, uprop->D[vp][nu][lu][cc3][cc2])
														, qcd_CONJ(qcd_CMUL(qcd_CSUB(PROJECTOR[0][ku][ku]
															, qcd_CMUL(imag, qcd_GAMMA[5][ku][ku])), usprop->D[v][ku][lu][cc1][cc2])))));

												block[v3] = qcd_CSUB(block[v3]
													, qcd_CMUL(wilsonline[v][cc1][cc3], qcd_CMUL(qcd_CMUL(cg5cg5m_val[ctrm2]
														, dprop->D[vp][nu][lu][cc3][cc2])
														, qcd_CONJ(qcd_CMUL(qcd_CADD(PROJECTOR[0][ku][ku]
															, qcd_CMUL(imag, qcd_GAMMA[5][ku][ku])), dsprop->D[v][ku][lu][cc1][cc2])))));

											}
										}
									}//space loop
								}//color2 loop    
							}//color1 loop

						}
			}


			if (geo->myid == 0)
			{
				// fprintf(fp_corr_p, "%+i %+i %+i %+i %+i", wl, mom[j][0], mom[j][1], mom[j][2],wl);
				fprintf(fp_corr_p, "%d  %d ", wl, t);

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

			MPI_Reduce(&(corr.re), &(corr2.re), 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			if (geo->myid == 0)
			{
				fprintf(fp_corr_p, "%+e %+e\n", corr2.re*0.5, corr2.im*0.5);
			}



		}//end lt inside local block condition
	}
	free(block);
	free(wilsonline);


}

void threepy(qcd_propagator *uprop, qcd_propagator *dprop, qcd_propagator*usprop, qcd_propagator *dsprop, \
	qcd_gaugeField *ustout, qcd_geometry *geo, qcd_int_4 x_src[4], qcd_int_4 t_start, qcd_int_4 t_stop, \
	qcd_int_4 w_length, FILE *fp_corr_p, qcd_int_4 mom[3])
{

	qcd_int_4 mu, ku, nu, lu;
	qcd_int_4 x, y, z, t, lx, ly, lz, lt;
	qcd_int_4 cc1, cc2, cc3, wl;
	qcd_int_8 v, v3, vm, vp;
	qcd_int_8 ctrm, ctrm2;

	qcd_real_8 tmp;
	qcd_complex_16 C, C2;
	qcd_complex_16 corr, corr2;
	qcd_int_2 cg5cg5m_ind[16 * 16][4];
	qcd_complex_16 cg5cg5m_val[16 * 16];
	qcd_complex_16 imag = { 0,1 };


	qcd_complex_16 *block;
	qcd_complex_16(*wilsonline)[3][3];
	qcd_complex_16 tempp[3][3];


	block = (qcd_complex_16*)malloc(geo->lV3 * sizeof(qcd_complex_16));
	wilsonline = malloc(geo->lV * 9 * sizeof(qcd_complex_16));

	if (wilsonline == NULL)
	{
		fprintf(stderr, "process %i: Error in wilsonline, Out of memory\n", geo->myid); exit(EXIT_FAILURE);
	}
	else
	{
		if (geo->myid == 0) printf("the wilsonline init sucessfully\n");
	}


	ctrm = 0;
	for (mu = 0; mu < 4; mu++)
		for (ku = 0; ku < 4; ku++)
			for (nu = 0; nu < 4; nu++)
			{
				C = qcd_CMUL(qcd_GAMMA[0][mu][nu], qcd_GAMMA[5][mu][ku]);
				if (qcd_NORM(C) > 1e-3)
				{
					cg5cg5m_val[ctrm].re = C.re;
					cg5cg5m_val[ctrm].im = C.im;
					cg5cg5m_ind[ctrm][0] = mu;
					cg5cg5m_ind[ctrm][1] = nu;
					cg5cg5m_ind[ctrm][2] = ku;
					ctrm++;
				}
			}

	for (v3 = 0; v3 < geo->lV; v3++)
	{
		qcd_unit3x3(wilsonline[v3]);
	}

	for (wl = 0; wl <= w_length; wl++)
	{
		for (t = t_start; t <= t_stop; t++)
		{
			lt = ((t + x_src[0]) % geo->L[0]) - geo->Pos[0] * geo->lL[0];

			for (v3 = 0; v3 < geo->lV3; v3++)   //set blocks to zero
				block[v3] = (qcd_complex_16) { 0, 0 };

			if (geo->myid == 0) printf(" calculating wilesonline %i time %i\n", wl, t);

			if (lt >= 0 && lt < geo->lL[0])  //inside the local lattice, otherwise nothing to calculate
			{

				for (lx = 0; lx < geo->lL[1]; lx++)
					for (ly = 0; ly < geo->lL[2]; ly++)
						for (lz = 0; lz < geo->lL[3]; lz++)
						{
							vp = qcd_LEXIC(lt, lx, (ly + wl - 1) % geo->lL[2], lz, geo->lL);
							vm = qcd_LEXIC(lt, lx, (ly + wl) % geo->lL[2], lz, geo->lL);
							v = qcd_LEXIC(lt, lx, ly, lz, geo->lL);
							v3 = qcd_LEXIC0(lx, ly, lz, geo->lL);
							if (wl != 0)
							{
								qcd_mul3x3(tempp, wilsonline[v], ustout->D[vp][2]);
								qcd_copy3x3(wilsonline[v], tempp);
							}
							for (lu = 0; lu < 4; lu++)
							{
								for (cc1 = 0; cc1 < 3; cc1++)
								{
									for (cc2 = 0; cc2 < 3; cc2++)
									{
										for (cc3 = 0; cc3 < 3; cc3++)
										{
											for (ctrm2 = 0; ctrm2 < ctrm; ctrm2++)
											{
												mu = cg5cg5m_ind[ctrm2][0];
												nu = cg5cg5m_ind[ctrm2][1];
												ku = cg5cg5m_ind[ctrm2][2];




												block[v3] = qcd_CADD(block[v3]
													, qcd_CMUL(wilsonline[v][cc1][cc3], qcd_CMUL(qcd_CMUL(cg5cg5m_val[ctrm2]
														, uprop->D[vm][nu][lu][cc3][cc2])
														, qcd_CONJ(qcd_CMUL(qcd_CSUB(PROJECTOR[0][ku][ku]
															, qcd_CMUL(imag, qcd_GAMMA[5][ku][ku])), usprop->D[v][ku][lu][cc1][cc2])))));

												block[v3] = qcd_CSUB(block[v3]
													, qcd_CMUL(wilsonline[v][cc1][cc3], qcd_CMUL(qcd_CMUL(cg5cg5m_val[ctrm2]
														, dprop->D[vm][nu][lu][cc3][cc2])
														, qcd_CONJ(qcd_CMUL(qcd_CADD(PROJECTOR[0][ku][ku]
															, qcd_CMUL(imag, qcd_GAMMA[5][ku][ku])), dsprop->D[v][ku][lu][cc1][cc2])))));

											}
										}
									}//space loop
								}//color2 loop    
							}//color1 loop

						}
			}



			if (geo->myid == 0)
			{
				// fprintf(fp_corr_p, "%+i %+i %+i %+i %+i", wl, mom[j][0], mom[j][1], mom[j][2],wl);
				fprintf(fp_corr_p, "%d  %d ", wl, t);

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

			MPI_Reduce(&(corr.re), &(corr2.re), 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			if (geo->myid == 0)
			{
				fprintf(fp_corr_p, "%+e %+e\n", corr2.re*0.5, corr2.im*0.5);
			}



		}//end lt inside local block condition

	}//end t-loop   

	for (v3 = 0; v3 < geo->lV; v3++)
	{
		qcd_unit3x3(wilsonline[v3]);
	}

	for (wl = -1; wl >= -w_length; wl--)
	{
		for (t = t_start; t <= t_stop; t++)
		{
			lt = ((t + x_src[0]) % geo->L[0]) - geo->Pos[0] * geo->lL[0];

			for (v3 = 0; v3 < geo->lV3; v3++)   //set blocks to zero
				block[v3] = (qcd_complex_16) { 0, 0 };

			if (geo->myid == 0) printf("calculating wilesonline %i time %i\n", wl, t);

			if (lt >= 0 && lt < geo->lL[0])  //inside the local lattice, otherwise nothing to calculate
			{

				for (lx = 0; lx < geo->lL[1]; lx++)
					for (ly = 0; ly < geo->lL[2]; ly++)
						for (lz = 0; lz < geo->lL[3]; lz++)
						{
							vp = qcd_LEXIC(lt, lx, (ly + wl + geo->lL[2]) % geo->lL[2], lz, geo->lL);
							v = qcd_LEXIC(lt, lx, ly, lz, geo->lL);
							v3 = qcd_LEXIC0(lx, ly, lz, geo->lL);
							qcd_mulAdjoint3x3(tempp, wilsonline[v], ustout->D[vp][2]);
							qcd_copy3x3(wilsonline[v], tempp);

							for (lu = 0; lu < 4; lu++)
							{
								for (cc1 = 0; cc1 < 3; cc1++)
								{
									for (cc2 = 0; cc2 < 3; cc2++)
									{
										for (cc3 = 0; cc3 < 3; cc3++)
										{
											for (ctrm2 = 0; ctrm2 < ctrm; ctrm2++)
											{
												mu = cg5cg5m_ind[ctrm2][0];
												nu = cg5cg5m_ind[ctrm2][1];
												ku = cg5cg5m_ind[ctrm2][2];


												block[v3] = qcd_CADD(block[v3]
													, qcd_CMUL(wilsonline[v][cc1][cc3], qcd_CMUL(qcd_CMUL(cg5cg5m_val[ctrm2]
														, uprop->D[vp][nu][lu][cc3][cc2])
														, qcd_CONJ(qcd_CMUL(qcd_CSUB(PROJECTOR[0][ku][ku]
															, qcd_CMUL(imag, qcd_GAMMA[5][ku][ku])), usprop->D[v][ku][lu][cc1][cc2])))));

												block[v3] = qcd_CSUB(block[v3]
													, qcd_CMUL(wilsonline[v][cc1][cc3], qcd_CMUL(qcd_CMUL(cg5cg5m_val[ctrm2]
														, dprop->D[vp][nu][lu][cc3][cc2])
														, qcd_CONJ(qcd_CMUL(qcd_CADD(PROJECTOR[0][ku][ku]
															, qcd_CMUL(imag, qcd_GAMMA[5][ku][ku])), dsprop->D[v][ku][lu][cc1][cc2])))));

											}
										}
									}//space loop
								}//color2 loop    
							}//color1 loop

						}
			}


			if (geo->myid == 0)
			{
				// fprintf(fp_corr_p, "%+i %+i %+i %+i %+i", wl, mom[j][0], mom[j][1], mom[j][2],wl);
				fprintf(fp_corr_p, "%d  %d ", wl, t);

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

			MPI_Reduce(&(corr.re), &(corr2.re), 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			if (geo->myid == 0)
			{
				fprintf(fp_corr_p, "%+e %+e\n", corr2.re*0.5, corr2.im*0.5);
			}



		}//end lt inside local block condition
	}
	free(block);
	free(wilsonline);


}




void threepz(qcd_propagator *uprop, qcd_propagator *dprop, qcd_propagator*usprop, qcd_propagator *dsprop, \
	qcd_gaugeField *ustout, qcd_geometry *geo, qcd_int_4 x_src[4], qcd_int_4 t_start, qcd_int_4 t_stop, \
	qcd_int_4 w_length, FILE *fp_corr_p, qcd_int_4 mom[3])
{

	qcd_int_4 mu, ku, nu, lu;
	qcd_int_4 x, y, z, t, lx, ly, lz, lt;
	qcd_int_4 cc1, cc2, cc3, wl;
	qcd_int_8 v, v3, vm, vp;
	qcd_int_8 ctrm, ctrm2;

	qcd_real_8 tmp;
	qcd_complex_16 C, C2;
	qcd_complex_16 corr, corr2;
	qcd_int_2 cg5cg5m_ind[16 * 16][4];
	qcd_complex_16 cg5cg5m_val[16 * 16];
	qcd_complex_16 imag = { 0,1 };


	qcd_complex_16 *block;
	qcd_complex_16(*wilsonline)[3][3];
	qcd_complex_16 tempp[3][3];


	block = (qcd_complex_16*)malloc(geo->lV3 * sizeof(qcd_complex_16));
	wilsonline = malloc(geo->lV * 9 * sizeof(qcd_complex_16));

	if (wilsonline == NULL)
	{
		fprintf(stderr, "process %i: Error in wilsonline, Out of memory\n", geo->myid); exit(EXIT_FAILURE);
	}
	else
	{
		if (geo->myid == 0) printf("the wilsonline init sucessfully\n");
	}


	ctrm = 0;
	for (mu = 0; mu < 4; mu++)
		for (ku = 0; ku < 4; ku++)
			for (nu = 0; nu < 4; nu++)
			{
				C = qcd_CMUL(qcd_GAMMA[0][mu][nu], qcd_GAMMA[5][mu][ku]);
				if (qcd_NORM(C) > 1e-3)
				{
					cg5cg5m_val[ctrm].re = C.re;
					cg5cg5m_val[ctrm].im = C.im;
					cg5cg5m_ind[ctrm][0] = mu;
					cg5cg5m_ind[ctrm][1] = nu;
					cg5cg5m_ind[ctrm][2] = ku;
					ctrm++;
				}
			}

	for (v3 = 0; v3 < geo->lV; v3++)
	{
		qcd_unit3x3(wilsonline[v3]);
	}

	for (wl = 0; wl <= w_length; wl++)
	{
		for (t = t_start; t <= t_stop; t++)
		{
			lt = ((t + x_src[0]) % geo->L[0]) - geo->Pos[0] * geo->lL[0];

			for (v3 = 0; v3 < geo->lV3; v3++)   //set blocks to zero
				block[v3] = (qcd_complex_16) { 0, 0 };

			if (geo->myid == 0) printf(" calculating wilesonline %i time %i\n", wl, t);

			if (lt >= 0 && lt < geo->lL[0])  //inside the local lattice, otherwise nothing to calculate
			{

				for (lx = 0; lx < geo->lL[1]; lx++)
					for (ly = 0; ly < geo->lL[2]; ly++)
						for (lz = 0; lz < geo->lL[3]; lz++)
						{
							vp = qcd_LEXIC(lt, lx, ly, (lz + wl - 1) % geo->lL[3], geo->lL);
							vm = qcd_LEXIC(lt, lx, ly, (lz + wl) % geo->lL[3], geo->lL);
							v = qcd_LEXIC(lt, lx, ly, lz, geo->lL);
							v3 = qcd_LEXIC0(lx, ly, lz, geo->lL);
							if (wl != 0)
							{
								qcd_mul3x3(tempp, wilsonline[v], ustout->D[vp][3]);
								qcd_copy3x3(wilsonline[v], tempp);
							}
							for (lu = 0; lu < 4; lu++)
							{
								for (cc1 = 0; cc1 < 3; cc1++)
								{
									for (cc2 = 0; cc2 < 3; cc2++)
									{
										for (cc3 = 0; cc3 < 3; cc3++)
										{
											for (ctrm2 = 0; ctrm2 < ctrm; ctrm2++)
											{
												mu = cg5cg5m_ind[ctrm2][0];
												nu = cg5cg5m_ind[ctrm2][1];
												ku = cg5cg5m_ind[ctrm2][2];




												block[v3] = qcd_CADD(block[v3]
													, qcd_CMUL(wilsonline[v][cc1][cc3], qcd_CMUL(qcd_CMUL(cg5cg5m_val[ctrm2]
														, uprop->D[vm][nu][lu][cc3][cc2])
														, qcd_CONJ(qcd_CMUL(qcd_CSUB(PROJECTOR[0][ku][ku]
															, qcd_CMUL(imag, qcd_GAMMA[5][ku][ku])), usprop->D[v][ku][lu][cc1][cc2])))));

												block[v3] = qcd_CSUB(block[v3]
													, qcd_CMUL(wilsonline[v][cc1][cc3], qcd_CMUL(qcd_CMUL(cg5cg5m_val[ctrm2]
														, dprop->D[vm][nu][lu][cc3][cc2])
														, qcd_CONJ(qcd_CMUL(qcd_CADD(PROJECTOR[0][ku][ku]
															, qcd_CMUL(imag, qcd_GAMMA[5][ku][ku])), dsprop->D[v][ku][lu][cc1][cc2])))));

											}
										}
									}//space loop
								}//color2 loop    
							}//color1 loop

						}
			}



			if (geo->myid == 0)
			{
				// fprintf(fp_corr_p, "%+i %+i %+i %+i %+i", wl, mom[j][0], mom[j][1], mom[j][2],wl);
				fprintf(fp_corr_p, "%d  %d ", wl, t);

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

			MPI_Reduce(&(corr.re), &(corr2.re), 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			if (geo->myid == 0)
			{
				fprintf(fp_corr_p, "%+e %+e\n", corr2.re*0.5, corr2.im*0.5);
			}



		}//end lt inside local block condition

	}//end t-loop   

	for (v3 = 0; v3 < geo->lV; v3++)
	{
		qcd_unit3x3(wilsonline[v3]);
	}

	for (wl = -1; wl >= -w_length; wl--)
	{
		for (t = t_start; t <= t_stop; t++)
		{
			lt = ((t + x_src[0]) % geo->L[0]) - geo->Pos[0] * geo->lL[0];

			for (v3 = 0; v3 < geo->lV3; v3++)   //set blocks to zero
				block[v3] = (qcd_complex_16) { 0, 0 };

			if (geo->myid == 0) printf("calculating wilesonline %i time %i\n", wl, t);

			if (lt >= 0 && lt < geo->lL[0])  //inside the local lattice, otherwise nothing to calculate
			{

				for (lx = 0; lx < geo->lL[1]; lx++)
					for (ly = 0; ly < geo->lL[2]; ly++)
						for (lz = 0; lz < geo->lL[3]; lz++)
						{
							vp = qcd_LEXIC(lt, lx, ly, (lz + wl + geo->lL[3]) % geo->lL[3], geo->lL);
							v = qcd_LEXIC(lt, lx, ly, lz, geo->lL);
							v3 = qcd_LEXIC0(lx, ly, lz, geo->lL);
							qcd_mulAdjoint3x3(tempp, wilsonline[v], ustout->D[vp][3]);
							qcd_copy3x3(wilsonline[v], tempp);

							for (lu = 0; lu < 4; lu++)
							{
								for (cc1 = 0; cc1 < 3; cc1++)
								{
									for (cc2 = 0; cc2 < 3; cc2++)
									{
										for (cc3 = 0; cc3 < 3; cc3++)
										{
											for (ctrm2 = 0; ctrm2 < ctrm; ctrm2++)
											{
												mu = cg5cg5m_ind[ctrm2][0];
												nu = cg5cg5m_ind[ctrm2][1];
												ku = cg5cg5m_ind[ctrm2][2];


												block[v3] = qcd_CADD(block[v3]
													, qcd_CMUL(wilsonline[v][cc1][cc3], qcd_CMUL(qcd_CMUL(cg5cg5m_val[ctrm2]
														, uprop->D[vp][nu][lu][cc3][cc2])
														, qcd_CONJ(qcd_CMUL(qcd_CSUB(PROJECTOR[0][ku][ku]
															, qcd_CMUL(imag, qcd_GAMMA[5][ku][ku])), usprop->D[v][ku][lu][cc1][cc2])))));

												block[v3] = qcd_CSUB(block[v3]
													, qcd_CMUL(wilsonline[v][cc1][cc3], qcd_CMUL(qcd_CMUL(cg5cg5m_val[ctrm2]
														, dprop->D[vp][nu][lu][cc3][cc2])
														, qcd_CONJ(qcd_CMUL(qcd_CADD(PROJECTOR[0][ku][ku]
															, qcd_CMUL(imag, qcd_GAMMA[5][ku][ku])), dsprop->D[v][ku][lu][cc1][cc2])))));

											}
										}
									}//space loop
								}//color2 loop    
							}//color1 loop

						}
			}


			if (geo->myid == 0)
			{
				// fprintf(fp_corr_p, "%+i %+i %+i %+i %+i", wl, mom[j][0], mom[j][1], mom[j][2],wl);
				fprintf(fp_corr_p, "%d  %d ", wl, t);

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

			MPI_Reduce(&(corr.re), &(corr2.re), 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			if (geo->myid == 0)
			{
				fprintf(fp_corr_p, "%+e %+e\n", corr2.re*0.5, corr2.im*0.5);
			}



		}//end lt inside local block condition
	}
	free(block);
	free(wilsonline);


}
