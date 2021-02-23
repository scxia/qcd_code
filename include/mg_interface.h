
#ifndef H_DDalphaAMG  
#define H_DDalphaAMG 1
  


#include <DDalphaAMG.h> 

extern qcd_uint_2 P[4];
extern qcd_uint_2 L[4];
extern qcd_geometry geo;
extern qcd_int_4 myid;

#define max(a,b) \
  ({ __typeof__ (a) _a = (a); \
  __typeof__ (b) _b = (b); \
  _a > _b ? _a : _b; })

static inline int conf_index_fct(int t, int z, int y, int x, int mu) { //T Z Y X
  int pos = qcd_LEXIC(t,x,y,z,geo.lL);
  int size_per_pos = 4*3*3*2; //link*SU(3)*complex
  int dir = (mu==1)?3:((mu==3)?1:mu); //swapping Z with X
  int size_per_dir = 3*3*2;   

  return pos*size_per_pos + dir*size_per_dir;
}

static inline int vector_index_fct(int t, int z, int y, int x) {
  int pos = qcd_LEXIC(t,x,y,z,geo.lL);
  int size_per_pos = 4*3*2; //link*SU(3)*complex
  
  return pos*size_per_pos;
}

static inline int Cart_rank(MPI_Comm comm, const int coords[], int *rank) {
  *rank = (int) qcd_LEXIC((P[0]+coords[0])%P[0],(P[1]+coords[3])%P[1],(P[2]+coords[2])%P[2],(P[3]+coords[1])%P[3],P);
  return *rank;
}

static inline int Cart_coords(MPI_Comm comm, int rank, int maxdims, int coords[]) {
  qcd_uint_2 tmp_p[4];
  qcd_antilexic(tmp_p, rank, P);
  coords[0]=tmp_p[0];
  coords[1]=tmp_p[3];
  coords[2]=tmp_p[2];
  coords[3]=tmp_p[1];
  return 1;
}

void mg_init(char* params, int params_len, double *theta) {
  
  DDalphaAMG_init mg_init;
  DDalphaAMG_parameters mg_params;
  DDalphaAMG_status status;
  mg_init.comm_cart = MPI_COMM_WORLD;
  mg_init.Cart_rank = Cart_rank;
  mg_init.Cart_coords = Cart_coords;
  mg_init.global_lattice[0] = L[0];
  mg_init.global_lattice[1] = L[3];
  mg_init.global_lattice[2] = L[2];
  mg_init.global_lattice[3] = L[1];
  mg_init.procs[0] = P[0];
  mg_init.procs[1] = P[3];
  mg_init.procs[2] = P[2];
  mg_init.procs[3] = P[1];
  mg_init.bc = 2; //(twisted)[if 2:multiply prop by exp(i*theta)]; 1 antiperiodic(not multiply)
  mg_init.theta[0] = theta[0]/M_PI;
  mg_init.theta[1] = theta[3]/M_PI;
  mg_init.theta[2] = theta[2]/M_PI;
  mg_init.theta[3] = theta[1]/M_PI;
  qcd_real_8 tol;
  
  /*sscanf(qcd_getParam("<tol>",params,params_len),"%lf", &tol);
    if(myid==0) printf("tolerance solver=%lf\n",tol); */

  sscanf(qcd_getParam("<block_txyz>",params,params_len),"%d %d %d %d",&(mg_init.block_lattice[0]), 
	  &(mg_init.block_lattice[3]), &(mg_init.block_lattice[2]), &(mg_init.block_lattice[1]));
#pragma omp parallel
   {
     mg_init.number_openmp_threads = omp_get_num_threads();
   }
  sscanf(qcd_getParam("<nlevels>",params,params_len),"%d", &(mg_init.number_of_levels));
  if(myid==0) printf("number of levels=%d\n",mg_init.number_of_levels);
  sscanf(qcd_getParam("<kappa>",params,params_len),"%lf",&(mg_init.kappa));
  if(myid==0) printf("kappa=%f\n",mg_init.kappa);
  sscanf(qcd_getParam("<mu>",params,params_len),"%lf",&(mg_init.mu));
  if(myid==0) printf("mu=%lf\n",mg_init.mu);
  sscanf(qcd_getParam("<csw>",params,params_len),"%lf",&(mg_init.csw));
  if(myid==0) printf("csw=%f\n",mg_init.csw);
  mg_init.init_file = NULL;
  mg_init.rnd_seeds = NULL;
  DDalphaAMG_initialize(&mg_init, &mg_params, &status);
 
  mg_params.setup_iterations[0] = 5;
  sscanf(qcd_getParam("<nvector>",params,params_len),"%d",&(mg_params.mg_basis_vectors[0]));
  mg_params.mg_basis_vectors[1] = max(28, mg_params.mg_basis_vectors[0]);
  for(int mu=0;mu<status.success; mu++)
    mg_params.mu_factor[mu]=1.;
  sscanf(qcd_getParam("<factor_cmu>",params,params_len),"%lf",&(mg_params.mu_factor[status.success-1]));
    
  mg_params.mu_odd_shift = 0;
  mg_params.mu_even_shift = 0;
  mg_params.mixed_precision = 1;
  mg_params.print = 3; //1 to have more printing
  mg_params.conf_index_fct = conf_index_fct;
  mg_params.vector_index_fct = vector_index_fct;
  DDalphaAMG_update_parameters(&mg_params, &status);

}

inline void mg_setup ( double *conf ) {

  DDalphaAMG_status status;

  DDalphaAMG_set_configuration( conf, &status );
  if(myid==0) printf("MG: plaquette = %e\n", status.info );
  
  if(myid==0) printf("MG: Running setup\n");
  DDalphaAMG_setup(&status);
  if(myid==0) printf("MG: Setup time %.2f sec (%.1f %% on coarse grid)\n", status.time,
		      100.*(status.coarse_time/status.time));

}

inline void mg_solve ( double *out, double *in, double tol ) {
  
  DDalphaAMG_status status;

  DDalphaAMG_solve( out, in, tol, &status );

  if(myid==0) printf("Solving time %.2f sec (%.1f %% on coarse grid)\n", status.time,
		     100.*(status.coarse_time/status.time));
  if(myid==0) printf("Total iterations on fine grid %d\n", status.iter_count);
  if(myid==0) printf("Total iterations on coarse grids %d\n", status.coarse_iter_count);
  if (!status.success && myid==0) printf("ERROR: not converged\n");
  if (!status.success) exit(0);
}	 

inline void mg_solve_up ( double *out, double *in, double tol ) {

  DDalphaAMG_parameters params;
  DDalphaAMG_status status;

  if(myid==0) printf("tolerance solver=%lf\n",tol);
  DDalphaAMG_get_parameters(&params);
  if( params.mu < 0 )
    DDalphaAMG_change_mu_sign(&status);

  mg_solve(out, in, tol);

}

inline void mg_solve_down ( double *out, double *in, double tol ) {

  DDalphaAMG_parameters params;
  DDalphaAMG_status status;

  DDalphaAMG_get_parameters(&params);
  if( params.mu > 0 )
    DDalphaAMG_change_mu_sign(&status);

  //TO DELETE
  DDalphaAMG_get_parameters(&params);
  if(myid==0) printf("MG: mu = %lf\n", params.mu);
  //END

  mg_solve(out, in, tol);

}

inline void mg_finalize () {
  DDalphaAMG_finalize();
}

#endif // !H_DDalphaAMG
