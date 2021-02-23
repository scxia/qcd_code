

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <qcd.h>
 

 

/* perform Ngauss iteration of gaussian smearing 
 * for a specific number of timeslices
 * parameter alpha.
 *
 * v-> (v+alpha*Hv) / (1+6*alpha)
 */
int qcd_gaussIteration3d_special_3pf(qcd_vector *v_u,qcd_vector *v_d,qcd_vector *v_s,qcd_vector *v_c,qcd_vector *v_su,qcd_vector *v_sd,
				     qcd_vector *v_ss,qcd_vector *v_sc,qcd_geometry *geo, qcd_gaugeField *u, qcd_real_8 alpha, qcd_uint_4 nsmear, qcd_uint_4 cur_time,
				     qcd_uint_4 number_tsinks, qcd_uint_4 after_tcur, qcd_uint_2 gauge_comm_flag)
{
   qcd_uint_8 i,j;
   qcd_uint_2 c1,mu,nu,b; 
   qcd_uint_4 x,y,z,b0,b1,b2,b3,tt=0;
   qcd_complex_16 tmp[3];
   qcd_real_8 nrm = 1.0/(1.0 + alpha*6.0);
   qcd_real_8 omega = alpha/(1.0 + alpha*6.0);
   qcd_real_8 *uu;
   qcd_real_8 *psi_fp[4],*psi_sp[4];
   qcd_real_8 upsi_fp[4][24],upsi_sp[4][24];
   qcd_real_8 udaggerpsi_fp[4][24],udaggerpsi_sp[4][24];
   qcd_real_8 *total_fp[4], *total_sp[4];
   qcd_uint_4 t;


   if(gauge_comm_flag == 1){
    qcd_communicateGaugeP(u);
    qcd_waitall(u->geo);
   }
    
   for(int ismear = 0 ; ismear < nsmear ; ismear++){

   qcd_communicateVectorPM(&v_u[0]);
      qcd_waitall(v_u[0].geo);
   qcd_communicateVectorPM(&v_d[0]);
      qcd_waitall(v_d[0].geo);
   qcd_communicateVectorPM(&v_s[0]);
      qcd_waitall(v_s[0].geo);
   qcd_communicateVectorPM(&v_c[0]);
      qcd_waitall(v_c[0].geo);
   qcd_communicateVectorPM(&v_su[0]); 
      qcd_waitall(v_su[0].geo);
   qcd_communicateVectorPM(&v_sd[0]); 
      qcd_waitall(v_sd[0].geo);
   qcd_communicateVectorPM(&v_ss[0]); 
      qcd_waitall(v_ss[0].geo);
   qcd_communicateVectorPM(&v_sc[0]);
      qcd_waitall(v_sc[0].geo);
   
//   if(geo->myid == 0){
//    printf("%+e %+e\n",v_u[0].D[geo->plus[10][2]][0][0].re,v_u[0].D[geo->plus[10][2]][0][0].im); 
//   }
      
  for(int it = 0; it < number_tsinks ; it++){
   t = ((it + cur_time + after_tcur)%geo->L[0]); 
   
   if( t>=geo->Pos[0]*geo->lL[0] && t<(geo->Pos[0]+1)*geo->lL[0])
   { 
      tt=t-geo->Pos[0]*geo->lL[0];
      for(z=0; z< geo->lL[3] ;z++)
      for(y=0; y< geo->lL[2] ;y++)
      for(x=0; x< geo->lL[1] ;x++)
      {
         i=qcd_LEXIC(tt,x,y,z,geo->lL);
	 
	 for(int spin = 0 ; spin < 4 ; spin++)
	   for(int color = 0 ; color < 3 ; color++){
	    v_u[1].D[i][spin][color] = (qcd_complex_16) {0,0};
	    v_d[1].D[i][spin][color] = (qcd_complex_16) {0,0};
	    v_s[1].D[i][spin][color] = (qcd_complex_16) {0,0};
	    v_c[1].D[i][spin][color] = (qcd_complex_16) {0,0};
	    v_su[1].D[i][spin][color] = (qcd_complex_16) {0,0};
	    v_sd[1].D[i][spin][color] = (qcd_complex_16) {0,0};
	    v_ss[1].D[i][spin][color] = (qcd_complex_16) {0,0};
	    v_sc[1].D[i][spin][color] = (qcd_complex_16) {0,0};
	   }
	   
         for(nu=1;nu<4;nu++)
         {
             uu  = (qcd_real_8*) &(u->D[i][nu][0][0].re);
             psi_fp[0] = (qcd_real_8*) &(v_u[0].D[geo->plus[i][nu]][0][0].re);
	     psi_fp[1] = (qcd_real_8*) &(v_d[0].D[geo->plus[i][nu]][0][0].re);
	     psi_fp[2] = (qcd_real_8*) &(v_s[0].D[geo->plus[i][nu]][0][0].re);
	     psi_fp[3] = (qcd_real_8*) &(v_c[0].D[geo->plus[i][nu]][0][0].re);
	     
	     psi_sp[0] = (qcd_real_8*) &(v_su[0].D[geo->plus[i][nu]][0][0].re);
	     psi_sp[1] = (qcd_real_8*) &(v_sd[0].D[geo->plus[i][nu]][0][0].re);
	     psi_sp[2] = (qcd_real_8*) &(v_ss[0].D[geo->plus[i][nu]][0][0].re);
	     psi_sp[3] = (qcd_real_8*) &(v_sc[0].D[geo->plus[i][nu]][0][0].re);
	     
	     for(int ii=0 ; ii < 4 ; ii++){
             qcd_APPLY_U(uu,upsi_fp[ii],psi_fp[ii]);
             qcd_APPLY_U(uu,upsi_sp[ii],psi_sp[ii]);	     
	     }
	     

             uu  = (qcd_real_8*) &(u->D[geo->minus[i][nu]][nu][0][0].re);
	     psi_fp[0] = (qcd_real_8*) &(v_u[0].D[geo->minus[i][nu]][0][0].re);
	     psi_fp[1] = (qcd_real_8*) &(v_d[0].D[geo->minus[i][nu]][0][0].re);
	     psi_fp[2] = (qcd_real_8*) &(v_s[0].D[geo->minus[i][nu]][0][0].re);
	     psi_fp[3] = (qcd_real_8*) &(v_c[0].D[geo->minus[i][nu]][0][0].re);
	     
	     psi_sp[0] = (qcd_real_8*) &(v_su[0].D[geo->minus[i][nu]][0][0].re);
	     psi_sp[1] = (qcd_real_8*) &(v_sd[0].D[geo->minus[i][nu]][0][0].re);
	     psi_sp[2] = (qcd_real_8*) &(v_ss[0].D[geo->minus[i][nu]][0][0].re);
	     psi_sp[3] = (qcd_real_8*) &(v_sc[0].D[geo->minus[i][nu]][0][0].re);
	      	     
	     for(int ii=0 ; ii < 4 ; ii++){
             qcd_APPLY_U_DAGGER(uu,udaggerpsi_fp[ii],psi_fp[ii]);
             qcd_APPLY_U_DAGGER(uu,udaggerpsi_sp[ii],psi_sp[ii]);	     
	     }
	     
	     for(int ii=0 ; ii< 4 ; ii++)
	       for(int jj=0 ; jj< 24 ; jj++){
		upsi_fp[ii][jj] = omega * upsi_fp[ii][jj];
		upsi_sp[ii][jj] = omega * upsi_sp[ii][jj];
		udaggerpsi_fp[ii][jj] = omega * udaggerpsi_fp[ii][jj];
		udaggerpsi_sp[ii][jj] = omega * udaggerpsi_sp[ii][jj];
	       }

	     total_fp[0] = (qcd_real_8*) &(v_u[1].D[i][0][0].re);
	     total_fp[1] = (qcd_real_8*) &(v_d[1].D[i][0][0].re);
	     total_fp[2] = (qcd_real_8*) &(v_s[1].D[i][0][0].re);
	     total_fp[3] = (qcd_real_8*) &(v_c[1].D[i][0][0].re);
	     
	     total_sp[0] = (qcd_real_8*) &(v_su[1].D[i][0][0].re);
	     total_sp[1] = (qcd_real_8*) &(v_sd[1].D[i][0][0].re);
	     total_sp[2] = (qcd_real_8*) &(v_ss[1].D[i][0][0].re);
	     total_sp[3] = (qcd_real_8*) &(v_sc[1].D[i][0][0].re);
	     
	     for(int i=0 ; i < 4 ; i++){
	      qcd_SUM_UP_HOPP(total_fp[i],upsi_fp[i],udaggerpsi_fp[i]);
	      qcd_SUM_UP_HOPP(total_sp[i],upsi_sp[i],udaggerpsi_sp[i]);
	     }
	     
         } // end directions
         	   
      }//end loop
   }//end condition
   
      for(z=0; z< geo->lL[3] ;z++)
      for(y=0; y< geo->lL[2] ;y++)
      for(x=0; x< geo->lL[1] ;x++)
      {
         i=qcd_LEXIC(tt,x,y,z,geo->lL);
           for(int spin=0 ; spin < 4 ; spin++)
	   for(int color = 0; color < 3 ; color++){
             v_u[0].D[i][spin][color].re = nrm * v_u[0].D[i][spin][color].re + v_u[1].D[i][spin][color].re ;
	     v_u[0].D[i][spin][color].im = nrm * v_u[0].D[i][spin][color].im + v_u[1].D[i][spin][color].im ;
	     
	     v_d[0].D[i][spin][color].re = nrm * v_d[0].D[i][spin][color].re + v_d[1].D[i][spin][color].re ;
	     v_d[0].D[i][spin][color].im = nrm * v_d[0].D[i][spin][color].im + v_d[1].D[i][spin][color].im ;
	     
	     v_s[0].D[i][spin][color].re = nrm * v_s[0].D[i][spin][color].re + v_s[1].D[i][spin][color].re ;
	     v_s[0].D[i][spin][color].im = nrm * v_s[0].D[i][spin][color].im + v_s[1].D[i][spin][color].im ;
	     
	     v_c[0].D[i][spin][color].re = nrm * v_c[0].D[i][spin][color].re + v_c[1].D[i][spin][color].re ;
	     v_c[0].D[i][spin][color].im = nrm * v_c[0].D[i][spin][color].im + v_c[1].D[i][spin][color].im ;
	     
	     v_su[0].D[i][spin][color].re = nrm * v_su[0].D[i][spin][color].re + v_su[1].D[i][spin][color].re ;
	     v_su[0].D[i][spin][color].im = nrm * v_su[0].D[i][spin][color].im + v_su[1].D[i][spin][color].im ;
	     
	     v_sd[0].D[i][spin][color].re = nrm * v_sd[0].D[i][spin][color].re + v_sd[1].D[i][spin][color].re ;
	     v_sd[0].D[i][spin][color].im = nrm * v_sd[0].D[i][spin][color].im + v_sd[1].D[i][spin][color].im ;
	     
	     v_ss[0].D[i][spin][color].re = nrm * v_ss[0].D[i][spin][color].re + v_ss[1].D[i][spin][color].re ;
	     v_ss[0].D[i][spin][color].im = nrm * v_ss[0].D[i][spin][color].im + v_ss[1].D[i][spin][color].im ;
	     
	     v_sc[0].D[i][spin][color].re = nrm * v_sc[0].D[i][spin][color].re + v_sc[1].D[i][spin][color].re ;
	     v_sc[0].D[i][spin][color].im = nrm * v_sc[0].D[i][spin][color].im + v_sc[1].D[i][spin][color].im ;
	   }
      }
   
  } // close time
  
 } // close smearing loop
   
   return 0;
} 

int qcd_gaussIteration3d_special_2pf(qcd_vector *v_u,qcd_vector *v_d,qcd_vector *v_s,qcd_vector *v_c,
                                     qcd_geometry *geo, qcd_gaugeField *u, qcd_real_8 alpha, qcd_uint_4 nsmear, qcd_uint_4 cur_time,
                                     qcd_uint_4 number_tsinks, qcd_uint_4 after_tcur, qcd_uint_2 gauge_comm_flag)
{
   qcd_uint_8 i,j;
   qcd_uint_2 c1,mu,nu,b; 
   qcd_uint_4 x,y,z,b0,b1,b2,b3,tt=0;
   qcd_complex_16 tmp[3];
   qcd_real_8 nrm = 1.0/(1.0 + alpha*6.0);
   qcd_real_8 omega = alpha/(1.0 + alpha*6.0);
   qcd_real_8 *uu;
   qcd_real_8 *psi_fp[4];
   qcd_real_8 upsi_fp[4][24];
   qcd_real_8 udaggerpsi_fp[4][24];
   qcd_real_8 *total_fp[4];
   qcd_uint_4 t;


   if(gauge_comm_flag == 1){
    qcd_communicateGaugeP(u);
    qcd_waitall(u->geo);
   }
    
   for(int ismear = 0 ; ismear < nsmear ; ismear++){

   qcd_communicateVectorPM(&v_u[0]);
      qcd_waitall(v_u[0].geo);
   qcd_communicateVectorPM(&v_d[0]);
      qcd_waitall(v_d[0].geo);
   qcd_communicateVectorPM(&v_s[0]);
      qcd_waitall(v_s[0].geo);
   qcd_communicateVectorPM(&v_c[0]);
      qcd_waitall(v_c[0].geo);
   
//   if(geo->myid == 0){
//    printf("%+e %+e\n",v_u[0].D[geo->plus[10][2]][0][0].re,v_u[0].D[geo->plus[10][2]][0][0].im); 
//   }
      
  for(int it = 0; it < number_tsinks ; it++){
   t = ((it + cur_time + after_tcur)%geo->L[0]); 
   
   if( t>=geo->Pos[0]*geo->lL[0] && t<(geo->Pos[0]+1)*geo->lL[0])
   { 
      tt=t-geo->Pos[0]*geo->lL[0];
      for(z=0; z< geo->lL[3] ;z++)
      for(y=0; y< geo->lL[2] ;y++)
      for(x=0; x< geo->lL[1] ;x++)
      {
         i=qcd_LEXIC(tt,x,y,z,geo->lL);
         
         for(int spin = 0 ; spin < 4 ; spin++)
           for(int color = 0 ; color < 3 ; color++){
            v_u[1].D[i][spin][color] = (qcd_complex_16) {0,0};
            v_d[1].D[i][spin][color] = (qcd_complex_16) {0,0};
            v_s[1].D[i][spin][color] = (qcd_complex_16) {0,0};
            v_c[1].D[i][spin][color] = (qcd_complex_16) {0,0};
           }
           
         for(nu=1;nu<4;nu++)
         {
             uu  = (qcd_real_8*) &(u->D[i][nu][0][0].re);
             psi_fp[0] = (qcd_real_8*) &(v_u[0].D[geo->plus[i][nu]][0][0].re);
             psi_fp[1] = (qcd_real_8*) &(v_d[0].D[geo->plus[i][nu]][0][0].re);
             psi_fp[2] = (qcd_real_8*) &(v_s[0].D[geo->plus[i][nu]][0][0].re);
             psi_fp[3] = (qcd_real_8*) &(v_c[0].D[geo->plus[i][nu]][0][0].re);
             
             
             for(int ii=0 ; ii < 4 ; ii++){
             qcd_APPLY_U(uu,upsi_fp[ii],psi_fp[ii]);
             }
             

             uu  = (qcd_real_8*) &(u->D[geo->minus[i][nu]][nu][0][0].re);
             psi_fp[0] = (qcd_real_8*) &(v_u[0].D[geo->minus[i][nu]][0][0].re);
             psi_fp[1] = (qcd_real_8*) &(v_d[0].D[geo->minus[i][nu]][0][0].re);
             psi_fp[2] = (qcd_real_8*) &(v_s[0].D[geo->minus[i][nu]][0][0].re);
             psi_fp[3] = (qcd_real_8*) &(v_c[0].D[geo->minus[i][nu]][0][0].re);
             
                     
             for(int ii=0 ; ii < 4 ; ii++){
             qcd_APPLY_U_DAGGER(uu,udaggerpsi_fp[ii],psi_fp[ii]);
             }
             
             for(int ii=0 ; ii< 4 ; ii++)
               for(int jj=0 ; jj< 24 ; jj++){
                upsi_fp[ii][jj] = omega * upsi_fp[ii][jj];
                udaggerpsi_fp[ii][jj] = omega * udaggerpsi_fp[ii][jj];
               }

             total_fp[0] = (qcd_real_8*) &(v_u[1].D[i][0][0].re);
             total_fp[1] = (qcd_real_8*) &(v_d[1].D[i][0][0].re);
             total_fp[2] = (qcd_real_8*) &(v_s[1].D[i][0][0].re);
             total_fp[3] = (qcd_real_8*) &(v_c[1].D[i][0][0].re);

             
             for(int i=0 ; i < 4 ; i++){
              qcd_SUM_UP_HOPP(total_fp[i],upsi_fp[i],udaggerpsi_fp[i]);
             }
             
         } // end directions
                   
      }//end loop
   }//end condition
   
      for(z=0; z< geo->lL[3] ;z++)
      for(y=0; y< geo->lL[2] ;y++)
      for(x=0; x< geo->lL[1] ;x++)
      {
         i=qcd_LEXIC(tt,x,y,z,geo->lL);
           for(int spin=0 ; spin < 4 ; spin++)
           for(int color = 0; color < 3 ; color++){
             v_u[0].D[i][spin][color].re = nrm * v_u[0].D[i][spin][color].re + v_u[1].D[i][spin][color].re ;
             v_u[0].D[i][spin][color].im = nrm * v_u[0].D[i][spin][color].im + v_u[1].D[i][spin][color].im ;
             
             v_d[0].D[i][spin][color].re = nrm * v_d[0].D[i][spin][color].re + v_d[1].D[i][spin][color].re ;
             v_d[0].D[i][spin][color].im = nrm * v_d[0].D[i][spin][color].im + v_d[1].D[i][spin][color].im ;
             
             v_s[0].D[i][spin][color].re = nrm * v_s[0].D[i][spin][color].re + v_s[1].D[i][spin][color].re ;
             v_s[0].D[i][spin][color].im = nrm * v_s[0].D[i][spin][color].im + v_s[1].D[i][spin][color].im ;
             
             v_c[0].D[i][spin][color].re = nrm * v_c[0].D[i][spin][color].re + v_c[1].D[i][spin][color].re ;
             v_c[0].D[i][spin][color].im = nrm * v_c[0].D[i][spin][color].im + v_c[1].D[i][spin][color].im ;
             
           }
      }
   
  } // close time
  
 } // close smearing loop
   
   return 0;
} 
