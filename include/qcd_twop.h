


#ifndef H_QCD_TWOP
#define H_QCD_TWOP 1

#include <stdio.h>

void twop_pi0_disconnect(qcd_vector * uprop, qcd_vector * dprop, qcd_geometry *geo, qcd_vector *vecsou, qcd_vector *vecsou2, \
	qcd_int_4 x_src[4], qcd_int_2 t_start, qcd_int_2 t_stop, qcd_int_4 mom[3], FILE *fp_corr_p);

void twop_delta(qcd_propagator * uprop, qcd_propagator* dprop, qcd_geometry *geo, qcd_int_4 x_src[4], qcd_int_2 t_start, \
	qcd_int_2 t_stop, qcd_int_4 mom[3], FILE *fp_corr_p);

void twop_proton(qcd_propagator * uprop, qcd_propagator* dprop, qcd_geometry *geo, qcd_int_4 x_src[4], qcd_int_2 t_start, \
	qcd_int_2 t_stop, qcd_int_4 mom[3], FILE *fp_corr_p);



#endif // !H_QCD_TWOP
