

#ifndef H_QCD_THREEP
#define H_QCD_THREEP 1 
#include <stdio.h>

void axial_charge(qcd_propagator *uprop, qcd_propagator *dprop, qcd_propagator*usprop, qcd_propagator *dsprop, \
	qcd_geometry *geo, qcd_int_4 x_src[4], qcd_int_4 t_start, qcd_int_4 t_stop, \
	FILE *fp_corr_p, qcd_int_4 mom[3]);




void threepx(qcd_propagator *uprop, qcd_propagator *dprop, qcd_propagator*usprop, qcd_propagator *dsprop, \
	qcd_gaugeField *ustout, qcd_geometry *geo, qcd_int_4 x_src[4], qcd_int_4 t_start, qcd_int_4 t_stop, \
	qcd_int_4 w_length, FILE *fp_corr_p, qcd_int_4 mom[3]);

void threepy(qcd_propagator *uprop, qcd_propagator *dprop, qcd_propagator*usprop, qcd_propagator *dsprop, \
	qcd_gaugeField *ustout, qcd_geometry *geo, qcd_int_4 x_src[4], qcd_int_4 t_start, qcd_int_4 t_stop, \
	qcd_int_4 w_length, FILE *fp_corr_p, qcd_int_4 mom[3]);


void threepz(qcd_propagator *uprop, qcd_propagator *dprop, qcd_propagator*usprop, qcd_propagator *dsprop, \
	qcd_gaugeField *ustout, qcd_geometry *geo, qcd_int_4 x_src[4], qcd_int_4 t_start, qcd_int_4 t_stop, \
	qcd_int_4 w_length, FILE *fp_corr_p, qcd_int_4 mom[3]);



#endif // !H_QCD_THREEP
