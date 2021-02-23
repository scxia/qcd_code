#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <math.h>
#include <qcd.h>
#include <mg_interface.h>


qcd_uint_2 P[4];
qcd_uint_2 L[4];
qcd_geometry geo;
qcd_int_4 myid;

void main(int argc, char* argv[])
{
	char   param_name[qcd_MAX_STRING_LENGTH];
	char   gauge_name[qcd_MAX_STRING_LENGTH];
	char   corr_p_name[qcd_MAX_STRING_LENGTH];
	char   three_p_name[qcd_MAX_STRING_LENGTH];
	char*  params;
	FILE *fp_corr_p,*fp_three_p;

	qcd_int_4 numprocs, params_len;
	qcd_int_4 i,k,is,lt,t,direction;


	//qcd_geometry geo;

	qcd_vector vecsou,vecout;
	qcd_propagator uprop, dprop;
	qcd_propagator upropsmear, dpropsmear;
	qcd_propagator suprop, sdprop;
	qcd_gaugeField u, uAPE, utemp,ustout;
	qcd_gaugeField* u_ptr, *uAPE_ptr, *utmp_ptr,*ustout_ptr;

	//qcd_uint_2 P[4];
	//qcd_uint_2 L[4];
	qcd_int_4 x_src[4];
	qcd_int_4 mom[3],momtf[3];
	qcd_int_4 t_start, t_stop,tsink,t_start3p,t_stop3p;
	qcd_int_4 w_length;
	qcd_real_8 theta[4] = { M_PI,0.,0.,0. };

	qcd_real_8  alpha, alphaAPE,keci,rhostout;
	qcd_int_4 nsmear, nsmearAPE,nsmearstout,stepstout;

	qcd_complex_16 phase_factor;
	qcd_real_8 accur;
	qcd_real_8 plaq;

	





	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);         // num. of processes taking part in the calculation
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);             // each process gets its ID
	//MPI_Get_processor_name(processor_name, &namelen);



	if (argc != 2)
	{
		if (myid == 0) fprintf(stderr, "No input file specified\n");
		exit(EXIT_FAILURE);
	}

	strcpy(param_name, argv[1]);


	if (myid == 0)
	{
		i = 0;
		printf("opening input file %s\n", param_name);
		params = qcd_getParams(param_name, &params_len);
		if (params == NULL)
		{
			i = 1;
		}
	}
	MPI_Bcast(&i, 1, MPI_INT, 0, MPI_COMM_WORLD);  //广播文件是否成功写入内存
	if (i == 1) exit(EXIT_FAILURE);
	MPI_Bcast(&params_len, 1, MPI_INT, 0, MPI_COMM_WORLD);   //广播内存大小
	if (myid != 0) params = (char*)malloc(params_len * sizeof(char));  //其他进程申请相同大小的内存
	MPI_Bcast(params, params_len, MPI_CHAR, 0, MPI_COMM_WORLD);  //广播数据




	sscanf(qcd_getParam("<processors_txyz>", params, params_len), "%hd %hd %hd %hd", &P[0], &P[1], &P[2], &P[3]);
	sscanf(qcd_getParam("<lattice_txyz>", params, params_len), "%hd %hd %hd %hd", &L[0], &L[1], &L[2], &L[3]);
	if (qcd_initGeometry(&geo, L, P, theta, myid, numprocs)) exit(EXIT_FAILURE);

	if (myid == 0) printf("Local lattice: %i x %i x %i x %i\n", geo.lL[0], geo.lL[1], geo.lL[2], geo.lL[3]);
	sscanf(qcd_getParam("<momentum_xyz>", params, params_len), "%d %d %d", &mom[0], &mom[1], &mom[2]);
	if (myid == 0) printf("Got momentum: %i %i %i \n", mom[0],mom[1],mom[2]);
	sscanf(qcd_getParam("<momentumtranf_xyz>", params, params_len), "%d %d %d", &momtf[0], &momtf[1], &momtf[2]);
	if (myid == 0) printf("Got momentum: %i %i %i \n", momtf[0], momtf[1], momtf[2]);
	sscanf(qcd_getParam("<alpha_gauss>", params, params_len), "%lf", &alpha);
	if (myid == 0) printf("Got alpha_gauss: %lf\n", alpha);  //高斯smearing的参数
	sscanf(qcd_getParam("<nsmear_gauss>", params, params_len), "%d", &nsmear);
	if (myid == 0) printf("Got nsmear_gauss: %d\n", nsmear);//高斯smearing的次数
	sscanf(qcd_getParam("<alpha_APE>", params, params_len), "%lf", &alphaAPE);
	if (myid == 0) printf("Got alpha_APE: %lf\n", alphaAPE);  //获取	APEsmearing的参数alpha
	sscanf(qcd_getParam("<nsmear_APE>", params, params_len), "%d", &nsmearAPE);
	if (myid == 0) printf("Got nsmear_APE: %d\n", nsmearAPE);       //获取APEsmearing的次数
	sscanf(qcd_getParam("<keci>", params, params_len), "%lf", &keci);
	if (myid == 0) printf("Got mommentum parameter keci: %lf\n", keci);  //高斯smearing的参数
	sscanf(qcd_getParam("<rho_stout>", params, params_len), "%lf", &rhostout);
	if (myid == 0) printf(" Got rho_stout: %lf\n", rhostout);
	sscanf(qcd_getParam("<nsmear_stout>", params, params_len), "%d", &nsmearstout);
	if (myid == 0) printf(" Got nsmear_stout: %d\n", nsmearstout);
	sscanf(qcd_getParam("<step_stout>", params, params_len), "%d", &stepstout);
	if (myid == 0) printf(" Got nsmear_stout: %d\n", stepstout);
	sscanf(qcd_getParam("<t>", params, params_len), "%d %d", &t_start, &t_stop);
	if (myid == 0) printf("Got sink time slices: %d ... %d\n", t_start, t_stop);
	sscanf(qcd_getParam("<t_3p>", params, params_len), "%d %d", &t_start3p, &t_stop3p);
	if (myid == 0) printf("Got sink time slices: %d ... %d\n", t_start3p, t_stop3p);
	sscanf(qcd_getParam("<tsink>", params, params_len), "%d ", &tsink);
	if (myid == 0) printf("Got sink time slices: %d \n", tsink);
	strcpy(gauge_name, qcd_getParam("<cfg_name>", params, params_len));
	if (myid == 0) printf("Got conf name: %s\n", gauge_name);
	strcpy(corr_p_name, qcd_getParam("<corr_name>", params, params_len));
	if (myid == 0) printf("Got twop output file name: %s\n", corr_p_name);
	strcpy(three_p_name, qcd_getParam("<three_name>", params, params_len));
	if (myid == 0) printf("Got threep output file name: %s\n", three_p_name);
	sscanf(qcd_getParam("<wilsonline>", params, params_len), "%d ", &w_length);
	if (myid == 0) printf("Got wilson length: %d \n", w_length);
	sscanf(qcd_getParam("<source_pos_txyz>", params, params_len), "%d %d %d %d", &x_src[0], &x_src[1], &x_src[2], &x_src[3]);
	if (myid == 0) printf("Got source coords: %d %d %d %d\n", x_src[0], x_src[1], x_src[2], x_src[3]);
	sscanf(qcd_getParam("<accur>", params, params_len), "%lf", &accur);
	if (myid == 0) printf("Got solver accurance: %+e\n", accur);
	sscanf(qcd_getParam("<direction>", params, params_len), "%d", &direction);
	if (myid == 0) printf("wilson direction: %d\n", direction);   // x->1 y->2 z->3 


	if (myid == 0)
	{
		k = 0;
		fp_corr_p = fopen(corr_p_name, "w");
		if (fp_corr_p == NULL)
		{
			printf("failed to open %s for writing twop\n", corr_p_name);
			k = 1;
		}
		fp_three_p = fopen(three_p_name, "w");
		if (fp_three_p == NULL)
		{
			printf("failed to open %s for writing threep \n", three_p_name);
			k = 1;
		}
	}
	MPI_Bcast(&k, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if (k == 1) exit(EXIT_FAILURE);



	qcd_initGaugeField(&u, &geo);
	qcd_initGaugeField(&uAPE, &geo);
	qcd_initGaugeField(&utemp, &geo);

	if (qcd_getGaugeField(gauge_name, qcd_GF_LIME, &u))//把组态gauge_name给u
	{
		fprintf(stderr, "process %i: Error reading gauge field!\n", myid);
		exit(EXIT_FAILURE);
	}
	if (myid == 0) printf("gauge-field loaded\n");
	plaq = qcd_calculatePlaquette(&u);
	if (myid == 0) printf("plaquette = %e\n", plaq);

	qcd_copyGaugeField(&utemp, &u);


	if (nsmear != 0)
	{
		u_ptr = &utemp;
		uAPE_ptr = &uAPE;
		for (i = 0; i < nsmearAPE; i++)    //做规范场的smearing
		{
			qcd_apeSmear3d(uAPE_ptr, u_ptr, alphaAPE);
			utmp_ptr = u_ptr; u_ptr = uAPE_ptr; uAPE_ptr = utmp_ptr;
		}
		utmp_ptr = u_ptr; u_ptr = uAPE_ptr; uAPE_ptr = utmp_ptr; //reverse the last swap. Also needed when nsmearAPE=0
		qcd_destroyGaugeField(u_ptr);
		uAPE = *uAPE_ptr;

		if (myid == 0) printf("gauge-field APE-smeared\n");
		plaq = qcd_calculatePlaquette(&uAPE);
		if (myid == 0) printf("plaquette = %e\n", plaq);
	}

	mg_init( params,  params_len, theta);
	mg_setup((double*) u.D);

	qcd_initVector(&vecsou, &geo);
	qcd_initVector(&vecout, &geo);
	qcd_initPropagator(&uprop, &geo);
	qcd_initPropagator(&dprop, &geo);
	qcd_initPropagator(&upropsmear, &geo);
	qcd_initPropagator(&dpropsmear, &geo);

	for (is = 0; is < 12; is++)
	{
		mksource(&vecsou, &uAPE, &geo, is, x_src, alpha, nsmear, keci, mom);

		mg_solve_up( (double*) vecout.D, (double*) vecsou.D, accur);
		qcd_copyPropagatorVector(&uprop, &vecout, is / 3, is % 3);
		for (i = 0; i < nsmear; i++)
		{
			if (qcd_momgaussIteration3dAll(&vecout, &uAPE, alpha, i == 0, mom, keci))
			{
				fprintf(stderr, "process %i: Error while smearing!\n", geo.myid);
				exit(EXIT_FAILURE);
			}
		}
		qcd_copyPropagatorVector(&upropsmear, &vecout, is / 3, is % 3);

		mg_solve_down((double *)vecout.D, (double *)vecsou.D, accur);
		qcd_copyPropagatorVector(&dprop, &vecout, is / 3, is % 3);
		for (i = 0; i < nsmear; i++)
		{
			if (qcd_momgaussIteration3dAll(&vecout, &uAPE, alpha, i == 0, mom, keci))
			{
				fprintf(stderr, "process %i: Error while smearing!\n", geo.myid);
				exit(EXIT_FAILURE);
			}
		}
		qcd_copyPropagatorVector(&dpropsmear, &vecout, is / 3, is % 3);

	}

	for (lt = 0; lt < geo.lL[0]; lt++)
	{
		t = (lt + geo.Pos[0] * geo.lL[0] - x_src[0] + geo.L[0]) % geo.L[0];    //t 是时间间隔
		phase_factor = (qcd_complex_16) { cos(theta[0] * t / geo.L[0]), sin(theta[0] * t / geo.L[0]) };

		qcd_mulPropagatorC3d(&uprop, phase_factor, lt);
		qcd_mulPropagatorC3d(&dprop, phase_factor, lt);
		qcd_mulPropagatorC3d(&upropsmear, phase_factor, lt);
		qcd_mulPropagatorC3d(&dpropsmear, phase_factor, lt);
	}

	qcd_tranformPropagatorPhysicalPlus(&uprop, &geo, 0, 0, geo.lL[0]);
	if (myid == 0) printf("up propagators transformed to physics basis \n");

	qcd_tranformPropagatorPhysicalMinus(&dprop, &geo, 0, 0, geo.lL[0]);
	if (myid == 0) printf("down propagators transformed to physics basis \n");

	qcd_tranformPropagatorPhysicalPlus(&upropsmear, &geo, 0, 0, geo.lL[0]);
	if (myid == 0) printf("up propagators transformed to physics basis \n");

	qcd_tranformPropagatorPhysicalMinus(&dpropsmear, &geo, 0, 0, geo.lL[0]);
	if (myid == 0) printf("down propagators transformed to physics basis \n");


	twop_proton(&upropsmear, &dpropsmear, &geo, x_src, t_start, t_stop, mom, fp_corr_p);

	if (myid == 0)
	{
		fclose(fp_corr_p);
	}

	qcd_initPropagator(&suprop, &geo);
	qcd_initPropagator(&sdprop, &geo);

	for (is = 0; is < 12; is++)
	{
		mksource_su_proton(&vecsou, &upropsmear, &dpropsmear, &uAPE, &geo, x_src, is, alpha, nsmear, keci, mom, tsink);
		mg_solve_down((double*)vecout.D, (double*)vecsou.D, accur);
		qcd_copyPropagatorVector(&suprop, &vecout, is / 3, is % 3);

		mksource_sd_proton(&vecsou, &upropsmear, &dpropsmear, &uAPE,& geo, x_src, is, alpha, nsmear, keci, mom, tsink);
		mg_solve_up((double*)vecout.D, (double*)vecsou.D, accur);
		qcd_copyPropagatorVector(&sdprop, &vecout, is / 3, is % 3);
		
	}

	qcd_destroyPropagator(&upropsmear);
	qcd_destroyPropagator(&dpropsmear);



	for (lt = 0; lt < geo.lL[0]; lt++)
	{
		t = (lt + geo.Pos[0] * geo.lL[0] - x_src[0] + geo.L[0]) % geo.L[0];    //t 是时间间隔
		phase_factor = (qcd_complex_16) { cos(theta[0] * (t-tsink) / geo.L[0]), sin(theta[0] * (t-tsink) / geo.L[0]) };

		qcd_mulPropagatorC3d(&suprop, phase_factor, lt);
		qcd_mulPropagatorC3d(&sdprop, phase_factor, lt);

	}


	qcd_initGaugeField(&ustout, &geo);
	u_ptr = &u;
	ustout_ptr = &ustout;



	for (i = 0; i < nsmearstout; i++)
	{
		if (i % stepstout == 0)
		{
			if (direction == 1)
			{
				threepx(&uprop, &dprop, &suprop, &sdprop, u_ptr, &geo, x_src, t_start3p, t_stop3p, w_length, fp_three_p, momtf);
			}
			if (direction == 2)
			{
				threepy(&uprop, &dprop, &suprop, &sdprop, u_ptr, &geo, x_src, t_start3p, t_stop3p, w_length, fp_three_p, momtf);
			}
			if (direction == 3)
			{
				threepz(&uprop, &dprop, &suprop, &sdprop, u_ptr, &geo, x_src, t_start3p, t_stop3p, w_length, fp_three_p, momtf);
			}

		}
		qcd_stoutsmearing(ustout_ptr, u_ptr, rhostout, 3);
		utmp_ptr = u_ptr; u_ptr = ustout_ptr; ustout_ptr = utmp_ptr;
	}

		
	


	if (myid == 0)
	{
		fclose(fp_three_p);
	}

	mg_finalize();
	qcd_destroyGaugeField(&u);
	qcd_destroyGaugeField(&uAPE);
	qcd_destroyGaugeField(&ustout);
	qcd_destroyVector(&vecsou);
	qcd_destroyVector(&vecout);
	qcd_destroyPropagator(&uprop);
	qcd_destroyPropagator(&dprop);
	qcd_destroyPropagator(&suprop);
	qcd_destroyPropagator(&sdprop);
	qcd_destroyGeometry(&geo);
	MPI_Finalize();


}