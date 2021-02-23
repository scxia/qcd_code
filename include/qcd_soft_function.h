
#ifndef H_QCD_SOFT_FUNCTION 
#define H_QCD_SOFT_FUNCTION 1

void twop_pion_soft_function(qcd_propagator *dprop_smear, qcd_propagator *upropsink_smear, qcd_geometry *geo, \
	qcd_int_4 x_src[4], qcd_int_4 x_sink[4], qcd_int_4 mom[3], FILE *fp_corr_p);

void form_factor_soft_function(qcd_propagator *dprop, qcd_propagator *upropsink, qcd_geometry *geo, \
	qcd_int_4 t_start, qcd_int_4 t_stop, qcd_int_4 x_src[4], qcd_int_4 x_sink[4], qcd_int_4 mom[3], qcd_int_4 momtf[3], FILE *fp_corr_p);

void pion_quasi_TMDWF_soft_function(qcd_propagator *dprop, qcd_geometry *geo, qcd_gaugeField *u, int x_src[4], FILE *fp_corr_p);

void pion_quasi_vacuum_value(qcd_geometry *geo, qcd_gaugeField *u, int x_src[4], FILE *fp_corr_p);
#endif