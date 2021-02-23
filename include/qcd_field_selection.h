#ifndef H_QCD_FIELDSELECTION
#define H_QCD_FIELDSELECTION 1


void twop_proton_field_selection(qcd_propagator * uprop, qcd_propagator* dprop, qcd_geometry *geo, qcd_int_4 x_src[4], qcd_int_2 t_start, \
	qcd_int_2 t_stop, qcd_int_4 mom[3], FILE *fp_corr_p, qcd_int_4(*random_list)[4], qcd_int_4 random_num);
void mksource_su_proton_field_selection(qcd_propagator *suprop, qcd_propagator *uprop, qcd_propagator *dprop, qcd_propagator *upropsink, qcd_propagator *dpropsink, qcd_gaugeField *uAPE, qcd_geometry *geo, \
	qcd_int_4 x_src[4], qcd_real_8 alpha, qcd_uint_4 nsmear, qcd_real_8 keci, qcd_int_4 mom[3], qcd_int_4 tsink, qcd_int_4 randompoint[4]);
void mksource_sd_proton_field_selection(qcd_propagator *sdprop, qcd_propagator *uprop, qcd_propagator *dprop, qcd_propagator *upropsink, qcd_propagator *dpropsink, qcd_gaugeField *uAPE, qcd_geometry *geo, \
	qcd_int_4 x_src[4], qcd_real_8 alpha, qcd_uint_4 nsmear, qcd_real_8 keci, qcd_int_4 mom[3], qcd_int_4 tsink, qcd_int_4 randompoint[4]);
void axial_charge_field_selection(qcd_propagator *uprop, qcd_propagator *dprop, qcd_propagator*usprop, qcd_propagator *dsprop, \
	qcd_geometry *geo, qcd_int_4 x_src[4], qcd_int_4 t_start, qcd_int_4 t_stop, FILE *fp_corr_p, qcd_int_4 mom[3], qcd_int_4(*random_list)[4], qcd_int_4 random_num, qcd_int_4 nsou, qcd_int_4 nsink);

#endif // !H_QCD_FIELDSELECTION




