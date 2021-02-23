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
	char   out_name_up_sp[qcd_MAX_STRING_LENGTH];
	char   out_name_up_ss[qcd_MAX_STRING_LENGTH];
	char   out_name_dn_sp[qcd_MAX_STRING_LENGTH];
	char   out_name_dn_ss[qcd_MAX_STRING_LENGTH];
	char   vec_names_up_sp[1024][qcd_MAX_STRING_LENGTH];
	char   vec_names_up_ss[1024][qcd_MAX_STRING_LENGTH];
	char   vec_names_dn_sp[1024][qcd_MAX_STRING_LENGTH];
	char   vec_names_dn_ss[1024][qcd_MAX_STRING_LENGTH];

	char*  params;
	FILE *pfile;

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

	
	int numb_sources=12;




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
	sscanf(qcd_getParam("<momentum_xyz>", params, params_len), "%d %d %d", &mom[0], &mom[1], &mom[2]);
	if (myid == 0) printf("Got momentum: %i %i %i \n", mom[0], mom[1], mom[2]);
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
	strcpy(gauge_name, qcd_getParam("<cfg_name>", params, params_len));
	if (myid == 0) printf("Got conf name: %s\n", gauge_name);
	sscanf(qcd_getParam("<source_pos_txyz>", params, params_len), "%d %d %d %d", &x_src[0], &x_src[1], &x_src[2], &x_src[3]);
	if (myid == 0) printf("Got source coords: %d %d %d %d\n", x_src[0], x_src[1], x_src[2], x_src[3]);
	sscanf(qcd_getParam("<accur>", params, params_len), "%lf", &accur);
	if (myid == 0) printf("Got solver accurance: %+e\n", accur);
	strcpy(out_name_up_sp, qcd_getParam("<propagator_up_sp>", params, params_len));    //source的输出文件列表
	if (myid == 0) printf("Got propagator name %s\n", out_name_up_sp);
	strcpy(out_name_up_ss, qcd_getParam("<propagator_up_ss>", params, params_len));    //source的输出文件列表
	if (myid == 0) printf("Got propagator name %s\n", out_name_up_ss);
	strcpy(out_name_dn_sp, qcd_getParam("<propagator_dn_sp>", params, params_len));    //source的输出文件列表
	if (myid == 0) printf("Got propagator name %s\n", out_name_dn_sp);
	strcpy(out_name_dn_ss, qcd_getParam("<propagator_dn_ss>", params, params_len));    //source的输出文件列表
	if (myid == 0) printf("Got propagator name %s\n", out_name_dn_ss);

//#######################################

	if ((pfile = fopen(out_name_up_sp, "r")) == NULL)
	{
		if (myid == 0) fprintf(stderr, "Error! Cannot open %s for reading.\n", out_name_up_sp);
		exit(EXIT_FAILURE);
	}
	for (i = 0; i < numb_sources; i++)
		fscanf(pfile, "%s\n", vec_names_up_sp[i]);
	fclose(pfile);

	if ((pfile = fopen(out_name_up_ss, "r")) == NULL)
	{
		if (myid == 0) fprintf(stderr, "Error! Cannot open %s for reading.\n", out_name_up_ss);
		exit(EXIT_FAILURE);
	}
	for (i = 0; i < numb_sources; i++)
		fscanf(pfile, "%s\n", vec_names_up_ss[i]);
	fclose(pfile);

	if ((pfile = fopen(out_name_dn_sp, "r")) == NULL)
	{
		if (myid == 0) fprintf(stderr, "Error! Cannot open %s for reading.\n", out_name_dn_sp);
		exit(EXIT_FAILURE);
	}
	for (i = 0; i < numb_sources; i++)
		fscanf(pfile, "%s\n", vec_names_dn_sp[i]);
	fclose(pfile);



	if ((pfile = fopen(out_name_dn_ss, "r")) == NULL)
	{
		if (myid == 0) fprintf(stderr, "Error! Cannot open %s for reading.\n", out_name_dn_ss);
		exit(EXIT_FAILURE);
	}
	for (i = 0; i < numb_sources; i++)
		fscanf(pfile, "%s\n", vec_names_dn_ss[i]);
	fclose(pfile);


	//##########################################


	qcd_initGaugeField(&u, &geo);
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
		qcd_initGaugeField(&uAPE, &geo);
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

		//qcd_initPropagator(&source,&geo);
	}


	mg_init( params,  params_len, theta);
	mg_setup((double*) u.D);

	qcd_initVector(&vecsou, &geo);
	qcd_initVector(&vecout, &geo);
	
	for (is = 0; is < 12; is++)
	{
		mksource(&vecsou, &uAPE, &geo, is, x_src, alpha, nsmear, keci, mom);

		mg_solve_up( (double*) vecout.D, (double*) vecsou.D, accur);
		qcd_writeVectorLime(vec_names_up_sp[is], qcd_PROP_LIME, &vecout);
		if (myid == 0) printf("up_sp propagator finished");
		for (i = 0; i < nsmear; i++)
		{
			if (qcd_momgaussIteration3dAll(&vecout, &uAPE, alpha, i == 0, mom, keci))
			{
				fprintf(stderr, "process %i: Error while smearing!\n", geo.myid);
				exit(EXIT_FAILURE);
			}
		}
		qcd_writeVectorLime(vec_names_up_ss[is], qcd_PROP_LIME, &vecout);
		if (myid == 0) printf("up_ss propagator finished");

		mg_solve_down((double *)vecout.D, (double *)vecsou.D, accur);
		qcd_writeVectorLime(vec_names_dn_sp[is], qcd_PROP_LIME, &vecout);
		if (myid == 0) printf("dn_sp propagator finished");

		for (i = 0; i < nsmear; i++)
		{
			if (qcd_momgaussIteration3dAll(&vecout, &uAPE, alpha, i == 0, mom, keci))
			{
				fprintf(stderr, "process %i: Error while smearing!\n", geo.myid);
				exit(EXIT_FAILURE);
			}
		}
		qcd_writeVectorLime(vec_names_dn_ss[is], qcd_PROP_LIME, &vecout);
		if (myid == 0) printf("dn_ss propagator finished");

		/*mksource(&vecsou, &u, &geo, is, x_src, 0, 0, keci, mom);

		mg_solve_up((double*)vecout.D, (double*)vecsou.D, accur);
		qcd_writeVectorLime(vec_names_up_pp[is], qcd_PROP_LIME, &vecout);
		if (myid == 0) printf("up_pp propagator finished");

		mg_solve_down((double *)vecout.D, (double *)vecsou.D, accur);
		qcd_writeVectorLime(vec_names_dn_pp[is], qcd_PROP_LIME, &vecout);
		if (myid == 0) printf("dn_pp propagator finished");*/
		


	}

	mg_finalize();
	qcd_destroyGaugeField(&u);
	qcd_destroyGaugeField(&uAPE);
	qcd_destroyVector(&vecsou);
	qcd_destroyVector(&vecout);
	qcd_destroyGeometry(&geo);
	MPI_Finalize();


}