
#ifndef H_QCD_SOURCE 
#define H_QCD_SOURCE 1

void mksource_su_delta(qcd_vector *souvec, qcd_propagator *uprop, qcd_propagator *dprop, qcd_gaugeField *uAPE, qcd_geometry *geo, \
	qcd_int_4 x_src[4], qcd_int_4 is, qcd_real_8 alpha, qcd_uint_4 nsmear, qcd_real_8 keci, qcd_int_4 mom[3], qcd_int_4 tsink);

void mksource_sd_delta(qcd_vector *souvec, qcd_propagator *uprop, qcd_propagator *dprop, qcd_gaugeField *uAPE, qcd_geometry *geo, \
	qcd_int_4 x_src[4], qcd_int_4 is, qcd_real_8 alpha, qcd_uint_4 nsmear, qcd_real_8 keci, qcd_int_4 mom[3], qcd_int_4 tsink);

void mksource_su_proton(qcd_vector *souvec, qcd_propagator *uprop, qcd_propagator *dprop, qcd_gaugeField *uAPE, qcd_geometry *geo, \
	qcd_int_4 x_src[4], qcd_int_4 is, qcd_real_8 alpha, qcd_uint_4 nsmear, qcd_real_8 keci, qcd_int_4 mom[3], qcd_int_4 tsink);

void mksource_sd_proton(qcd_vector *souvec, qcd_propagator *uprop, qcd_propagator *dprop, qcd_gaugeField *uAPE, qcd_geometry *geo, \
	qcd_int_4 x_src[4], qcd_int_4 is, qcd_real_8 alpha, qcd_uint_4 nsmear, qcd_real_8 keci, qcd_int_4 mom[3], qcd_int_4 tsink);

void mksource(qcd_vector *vec, qcd_gaugeField *uAPE, qcd_geometry *geo, qcd_int_4 is, qcd_int_4 x_src[4], qcd_real_8 alpha, qcd_uint_4 nsmear, qcd_real_8 keci, qcd_int_4 mom[3]);

void WallSource3D(qcd_vector* vec, qcd_int_4 mom[3], qcd_int_4 t, qcd_int_4 is);

void WallSource4D(qcd_vector* vec, qcd_int_4 mom[4], qcd_int_4 is);

void mksource_stochasticZ4(qcd_vector* vec, qcd_int_4 t, qcd_int_4 is);
void stochasticZ4(qcd_vector* des, qcd_vector* vec, qcd_int_4 mom[3], qcd_int_4 t, qcd_int_4 dis, qcd_int_4 vis);

#endif // !H_QCD_SOURCE 1
