#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <qcd.h>

void qcd_copyVectorPropagator_timerange(qcd_vector *vec, qcd_propagator *prop, qcd_uint_2 nu, qcd_uint_2 c2, qcd_uint_4 cur_time, qcd_uint_4 after_tcur, qcd_uint_4 number_tsinks)
{
  qcd_uint_8 i;
  qcd_uint_2 mu,c1;
  qcd_uint_4 t;
  qcd_uint_4 z,y,x;

  for(int it = 0; it < number_tsinks ; it++){
    t = ( (it + cur_time + after_tcur) % vec->geo->L[0]);

    for(z=0; z < vec->geo->lL[3] ;z++)
      for(y=0; y < vec->geo->lL[2] ;y++)
	for(x=0; x < vec->geo->lL[1] ;x++){
	  i=qcd_LEXIC(t,x,y,z,vec->geo->lL);
	  for(mu=0; mu<4; mu++)
	    for(c1=0; c1<3; c1++)
	      vec->D[i][mu][c1] = prop->D[i][mu][nu][c1][c2];
	}
  }
}

void qcd_copyPropagatorVector_timerange(qcd_propagator *prop, qcd_vector *vec, qcd_uint_2 nu, qcd_uint_2 c2, qcd_uint_4 cur_time, qcd_uint_4 after_tcur, qcd_uint_4 number_tsinks)
{
  qcd_uint_8 i;
  qcd_uint_2 mu,c1;
  qcd_uint_4 t;
  qcd_uint_4 z,y,x;

  for(int it = 0; it < number_tsinks ; it++){
    t = ( (it + cur_time + after_tcur) % vec->geo->L[0]);

    for(z=0; z < vec->geo->lL[3] ;z++)
      for(y=0; y < vec->geo->lL[2] ;y++)
	for(x=0; x < vec->geo->lL[1] ;x++){
          i=qcd_LEXIC(t,x,y,z,vec->geo->lL);
	  
	  for(mu=0; mu<4; mu++)
	    for(c1=0; c1<3; c1++)
	      prop->D[i][mu][nu][c1][c2] = vec->D[i][mu][c1];
	}
  }
}

void qcd_tranformPropagatorPhysicalPlus(qcd_propagator *prop, qcd_geometry *geo, qcd_uint_4 cur_time, qcd_uint_4 after_tcur, qcd_uint_4 number_tsinks){

  qcd_complex_16 onePg5[4][4];
  qcd_complex_16 one[4][4];
  qcd_complex_16 imag, C;
  qcd_int_2 onePg5Sq_ind[16*16][4];
  qcd_complex_16 onePg5Sq_val[16*16];
  qcd_uint_4 counter = 0;
  qcd_uint_2 t;
  qcd_uint_8 v;
  qcd_complex_16 tmp[4][4];
  qcd_uint_2 alpha, beta, gamma, delta;
 
  imag = (qcd_complex_16) {0,1};
 
  for(int i=0; i<4 ; i++)
    for(int j=0 ; j<4 ; j++)
      one[i][j] = (qcd_complex_16) {0,0};
  
  one[0][0] = (qcd_complex_16) {1,0};
  one[1][1] = (qcd_complex_16) {1,0};
  one[2][2] = (qcd_complex_16) {1,0};
  one[3][3] = (qcd_complex_16) {1,0};
 
  for(int i= 0 ; i < 4 ; i++)
    for(int j=0; j < 4 ; j++)
      onePg5[i][j] = qcd_CADD(one[i][j],qcd_CMUL(imag,qcd_GAMMA[5][i][j]));
  
  for(int alpha = 0 ; alpha < 4 ; alpha++)
    for(int beta = 0 ; beta < 4 ; beta++)
      for(int gamma = 0 ; gamma < 4 ; gamma++)
	for(int delta = 0 ; delta < 4 ; delta++)
	  {
	    C = qcd_CMUL(onePg5[alpha][beta], onePg5[gamma][delta]);
	    if(qcd_NORM(C)>1e-3){
	      onePg5Sq_val[counter].re = 0.5 * C.re;
	      onePg5Sq_val[counter].im = 0.5 * C.im;
	      onePg5Sq_ind[counter][0] = alpha ;
	      onePg5Sq_ind[counter][1] = beta ;
	      onePg5Sq_ind[counter][2] = gamma ;
	      onePg5Sq_ind[counter][3] = delta ;
	      counter++;
	    }
	  }
     
   
  for(int it = 0; it < number_tsinks ; it++){
    t = ( (it + cur_time + after_tcur) % geo->L[0]);

    for(int z=0; z < geo->lL[3] ;z++)
      for(int y=0; y < geo->lL[2] ;y++)
        for(int x=0; x < geo->lL[1] ;x++){
          v=qcd_LEXIC(t,x,y,z,geo->lL);
          
          for(int c1 = 0 ; c1 < 3 ; c1++)
	    for(int c2 = 0 ; c2 < 3 ; c2++){
            
	      for(int i = 0 ; i < 4 ; i++)
		for(int j = 0 ; j< 4 ; j++)
		  tmp[i][j] = (qcd_complex_16) {0,0};
            
	      for(int spins =0 ; spins < counter ; spins++){
		alpha = onePg5Sq_ind[spins][0];
		beta = onePg5Sq_ind[spins][1];
		gamma = onePg5Sq_ind[spins][2];
		delta = onePg5Sq_ind[spins][3];
              
		tmp[alpha][delta] = qcd_CADD(tmp[alpha][delta],qcd_CMUL(onePg5Sq_val[spins],prop->D[v][beta][gamma][c1][c2]));
	      } // close spins
            
	      for(int alpha = 0 ; alpha < 4 ; alpha++)
		for(int delta = 0 ; delta < 4 ; delta++)
		  prop->D[v][alpha][delta][c1][c2] = tmp[alpha][delta];
             
	    } // close colors
          
        } // close spatial local volume
        
  } // close time
  
 
}

void qcd_tranformPropagatorPhysicalMinus(qcd_propagator *prop, qcd_geometry *geo, qcd_uint_4 cur_time, qcd_uint_4 after_tcur, qcd_uint_4 number_tsinks){

  qcd_complex_16 oneMg5[4][4];
  qcd_complex_16 one[4][4];
  qcd_complex_16 imag, C;
  qcd_int_2 oneMg5Sq_ind[16*16][4];
  qcd_complex_16 oneMg5Sq_val[16*16];
  qcd_uint_4 counter = 0;
  qcd_uint_2 t;
  qcd_uint_8 v;
  qcd_complex_16 tmp[4][4];
  qcd_uint_2 alpha, beta, gamma, delta;
 
  imag = (qcd_complex_16) {0,1};
 
  for(int i=0; i<4 ; i++)
    for(int j=0 ; j<4 ; j++)
      one[i][j] = (qcd_complex_16) {0,0};
  
  one[0][0] = (qcd_complex_16) {1,0};
  one[1][1] = (qcd_complex_16) {1,0};
  one[2][2] = (qcd_complex_16) {1,0};
  one[3][3] = (qcd_complex_16) {1,0};
 
  for(int i= 0 ; i < 4 ; i++)
    for(int j=0; j < 4 ; j++)
      oneMg5[i][j] = qcd_CSUB(one[i][j],qcd_CMUL(imag,qcd_GAMMA[5][i][j]));
  
  for(int alpha = 0 ; alpha < 4 ; alpha++)
    for(int beta = 0 ; beta < 4 ; beta++)
      for(int gamma = 0 ; gamma < 4 ; gamma++)
	for(int delta = 0 ; delta < 4 ; delta++)
	  {
	    C = qcd_CMUL(oneMg5[alpha][beta], oneMg5[gamma][delta]);
	    if(qcd_NORM(C)>1e-3){
	      oneMg5Sq_val[counter].re = 0.5 * C.re;
	      oneMg5Sq_val[counter].im = 0.5 * C.im;
	      oneMg5Sq_ind[counter][0] = alpha ;
	      oneMg5Sq_ind[counter][1] = beta ;
	      oneMg5Sq_ind[counter][2] = gamma ;
	      oneMg5Sq_ind[counter][3] = delta ;
	      counter++;
	    }
	  }
     
   
  for(int it = 0; it < number_tsinks ; it++){
    t = ( (it + cur_time + after_tcur) % geo->L[0]);

    for(int z=0; z < geo->lL[3] ;z++)
      for(int y=0; y < geo->lL[2] ;y++)
        for(int x=0; x < geo->lL[1] ;x++){
          v=qcd_LEXIC(t,x,y,z,geo->lL);
          
          for(int c1 = 0 ; c1 < 3 ; c1++)
	    for(int c2 = 0 ; c2 < 3 ; c2++){
            
	      for(int i = 0 ; i < 4 ; i++)
		for(int j = 0 ; j< 4 ; j++)
		  tmp[i][j] = (qcd_complex_16) {0,0};
            
	      for(int spins =0 ; spins < counter ; spins++){
		alpha = oneMg5Sq_ind[spins][0];
		beta = oneMg5Sq_ind[spins][1];
		gamma = oneMg5Sq_ind[spins][2];
		delta = oneMg5Sq_ind[spins][3];
              
		tmp[alpha][delta] = qcd_CADD(tmp[alpha][delta],qcd_CMUL(oneMg5Sq_val[spins],prop->D[v][beta][gamma][c1][c2]));
	      } // close spins
            
	      for(int alpha = 0 ; alpha < 4 ; alpha++)
		for(int delta = 0 ; delta < 4 ; delta++)
		  prop->D[v][alpha][delta][c1][c2] = tmp[alpha][delta];
             
	    } // close colors
          
        } // close spatial local volume
        
  } // close time
  
 
}



char* startcode( qcd_int_4  *params_len,qcd_int_4 *myid,qcd_int_4 *numprocs,int argc, char* argv[])
{
	char   param_name[qcd_MAX_STRING_LENGTH];
	char* params;
	qcd_int_4 i;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, numprocs);         // num. of processes taking part in the calculation
	MPI_Comm_rank(MPI_COMM_WORLD, myid);             // each process gets its ID
	//MPI_Get_processor_name(processor_name, &namelen);
	if (argc != 2)
	{
		if (*myid == 0) fprintf(stderr, "No input file specified\n");
		exit(EXIT_FAILURE);
	}

	strcpy(param_name, argv[1]);


	if (*myid == 0)
	{
		i = 0;
		printf("opening input file %s\n", param_name);
		params = qcd_getParams(param_name, params_len);
		if (params == NULL)
		{
			i = 1;
		}
	}
	MPI_Bcast(&i, 1, MPI_INT, 0, MPI_COMM_WORLD);  //广播文件是否成功写入内存
	if (i == 1) exit(EXIT_FAILURE);
	MPI_Bcast(params_len, 1, MPI_INT, 0, MPI_COMM_WORLD);   //广播内存大小
	if (*myid != 0) params = (char*)malloc(*params_len * sizeof(char));  //其他进程申请相同大小的内存
	MPI_Bcast(params, *params_len, MPI_CHAR, 0, MPI_COMM_WORLD);  //广播数据
	return params;
}


void MomSmearing(qcd_propagator *prop,qcd_gaugeField *uAPE,qcd_int_4 mom[3],qcd_int_4 nsmear,qcd_real_8 alpha,qcd_real_8 keci)
{

	qcd_int_4 mu ,c1,i;
	qcd_vector vec;
	qcd_initVector(&vec,prop->geo);
	for (mu = 0; mu < 4; mu++)
		for (c1 = 0; c1 < 3; c1++)
		{
			qcd_copyVectorPropagator(&vec, prop, mu, c1);  
			for (i = 0; i < nsmear; i++)
			{
				if (qcd_momgaussIteration3dAll(&vec, uAPE, alpha, i == 0, mom, keci))
				{
					fprintf(stderr, "process %i: Error while smearing!\n", prop->geo->myid);
					exit(EXIT_FAILURE);
				}
			}
			qcd_copyPropagatorVector(prop, &vec, mu, c1);
			
		}
	if (prop->geo->myid == 0)
	{
		printf("MomSmearing finished \n");
	}
	qcd_destroyVector(&vec);
}

void APESmearing(qcd_gaugeField* u, qcd_int_4 nsmearAPE,qcd_real_8 alphaAPE)
{
	qcd_gaugeField uAPE;
	qcd_gaugeField* u_ptr, * uAPE_ptr, * utmp_ptr;
	qcd_int_4 i;
	qcd_real_8 plaq;

	qcd_initGaugeField(&uAPE, u->geo);
	uAPE_ptr = &uAPE;
	u_ptr = u;
	for (i = 0; i < nsmearAPE; i++)
	{
		qcd_apeSmear3d(uAPE_ptr, u_ptr, alphaAPE);
		utmp_ptr = u_ptr; u_ptr = uAPE_ptr; uAPE_ptr = utmp_ptr;
	}
	utmp_ptr = u_ptr; u_ptr = uAPE_ptr; uAPE_ptr = utmp_ptr; //reverse the last swap. Also needed when nsmearAPE=0
	qcd_destroyGaugeField(u_ptr);
	uAPE = *uAPE_ptr;

	u = &uAPE;

	plaq = qcd_calculatePlaquette(u);
	if (u->geo->myid == 0) printf("After APE smearing plaquette = %e\n", plaq);
}







