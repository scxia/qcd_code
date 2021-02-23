#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <math.h>
#include <qcd.h>
#include<qcd_soft_function.h>
#include<qcd_gaugefix.h>


qcd_uint_2 P[4];
qcd_uint_2 L[4];
qcd_geometry geo;
qcd_int_4 myid;

void main(int argc, char* argv[])
{

	char   param_name[qcd_MAX_STRING_LENGTH];
	char   gauge_name[qcd_MAX_STRING_LENGTH];
	char   gauge_name_gaugefix[qcd_MAX_STRING_LENGTH];
	char* params;
	qcd_int_4 numprocs, params_len;
	qcd_int_4 i;
	qcd_real_8 plaq;
	qcd_gaugeField u;
	qcd_gaugeField uc;

	qcd_real_8 theta[4] = { M_PI,0.,0.,0. };

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);         // num. of processes taking part in the calculation
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);


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

	strcpy(gauge_name, qcd_getParam("<cfg_name>", params, params_len));
	if (myid == 0) printf("Got conf name: %s\n", gauge_name);

	strcpy(gauge_name_gaugefix, qcd_getParam("<cfg_name_gaugefix>", params, params_len));
	if (myid == 0) printf("Got conf gauge fix name: %s\n", gauge_name_gaugefix);

	qcd_initGaugeField(&u, &geo);

	if (qcd_getGaugeField(gauge_name, qcd_GF_LIME, &u))
	{
		if (myid == 0) fprintf(stderr, "Error reading gauge field\n");
		exit(EXIT_FAILURE);
	}
	if (myid == 0) printf("gauge-field loaded\n");
	plaq = qcd_calculatePlaquette(&u);
	if (myid == 0) printf("plaquette = %e\n", plaq);


		
	gaugefix(Landau, &u, 1e-8, 5000, 1);


	plaq = qcd_calculatePlaquette(&u);
	if (myid == 0) printf("plaquette = %e\n", plaq);

	qcd_writeGaugeField(gauge_name_gaugefix, qcd_GF_LIME, &u);

}