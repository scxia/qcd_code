#ifndef H_QCD_SOFT_FUNCTION_WALLSOURCE 
#define H_QCD_SOFT_FUNCTION_WALLSOURCE 1
void PropagatorFourierTransform3D(qcd_complex_16 res2[4][4][3][3], qcd_propagator* u, qcd_int_4 mom[3], qcd_int_4 x_src[4], qcd_int_4 lt);
void PropagatorFourierTransform4D(qcd_complex_16 res2[4][4][3][3], qcd_propagator* u, qcd_int_4 mom[4], qcd_int_4 x_src[4]);
void BlockMatrixFourierTransform3D(qcd_complex_16 res2[4][4], qcd_geometry* geo, qcd_complex_16(*block)[4][4], \
	qcd_int_4 mom[3], qcd_int_4 x_src[4], qcd_int_4 lt);
void BlockFourierTransform3D(qcd_complex_16 *res2, qcd_geometry* geo, qcd_complex_16* block, qcd_int_4 mom[3], \
	qcd_int_4 x_src[4], qcd_int_4 lt);
void PionTwoPointWallSource(qcd_propagator* u1, qcd_propagator* u2, qcd_int_4 mom1[3], qcd_int_4 mom2[3], \
	qcd_int_4 tstart, qcd_int_4 tend, qcd_int_4 src_t, FILE* fp_corr_p);
void PionTwoPointMomSmearing(qcd_propagator* propp, qcd_propagator* propm, qcd_geometry* geo, qcd_int_4 tstart, qcd_int_4 tend, \
	qcd_int_4 mom[3], qcd_int_4 x_src[4], FILE* fp_corr_p);
void PionFormFactorWallSource(qcd_propagator *uprop1, qcd_propagator *uprop2, qcd_propagator *dprop1, qcd_propagator *dprop2, \
	qcd_geometry* geo, qcd_int_4 t_start, qcd_int_4 t_stop, qcd_int_4 momtf[3], qcd_int_4 tsource, FILE* fp_corr_p);
void StapleWilsonLine(qcd_complex_16 Staple[3][3], int lt, int lx, int ly, int lz, int l, int bprep, qcd_geometry* geo, qcd_gaugeField* u);
void PionQuasiWaveFunctionWallSource(qcd_propagator* prop1, qcd_propagator* prop2, qcd_gaugeField* u, qcd_geometry* geo,\
	qcd_int_4 mom[3], qcd_int_4 src_t, FILE* fp_corr_p);
void QuarkLocalOperatorRenormalization(qcd_propagator* prop, FILE* fp_corr_p);
void QuarkQuasiWaveFunctionWallSource(qcd_propagator* prop, qcd_gaugeField* u, qcd_geometry* geo, FILE* fp_corr_p);
void QuarkWaveFunctionFourierTransformation(qcd_propagator* prop, qcd_int_4 mom[4], FILE* fp_corr_p);
void PionTwoPointAxial(qcd_propagator* propp, qcd_propagator* propm, qcd_geometry* geo, qcd_int_4 tstart, qcd_int_4 tend, \
	qcd_int_4 mom[3], qcd_int_4 x_src[4], FILE* fp_corr_p);
#endif