#ifndef H_QCD_SMEARING_SPECIAL
#define H_QCD_SMEARING_SPECIAL 1
 
int qcd_gaussIteration3d_special_3pf(qcd_vector *v_u,qcd_vector *v_d,qcd_vector *v_s,qcd_vector *v_c,qcd_vector *v_su,
qcd_vector *v_sd, qcd_vector *v_ss,qcd_vector *v_sc,qcd_geometry *geo, qcd_gaugeField *u, qcd_real_8 alpha, qcd_uint_4 nsmear, 
				     qcd_uint_4 cur_time, qcd_uint_4 number_tsinks, qcd_uint_4 after_tcur, qcd_uint_2 gauge_comm_flag);

int qcd_gaussIteration3d_special_2pf(qcd_vector *v_u,qcd_vector *v_d,qcd_vector *v_s,qcd_vector *v_c
				     ,qcd_geometry *geo, qcd_gaugeField *u, qcd_real_8 alpha, qcd_uint_4 nsmear, 
				     qcd_uint_4 cur_time, qcd_uint_4 number_tsinks, qcd_uint_4 after_tcur, qcd_uint_2 gauge_comm_flag);


#endif
