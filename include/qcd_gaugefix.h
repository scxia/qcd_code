#ifndef H_QCD_GAUGEFIX
#define H_QCD_GAUGEFIX 1

#define Landau 0
#define Coulomb 1
double gauge_fix_action(qcd_gaugeField* u, qcd_int_4 OrT);
void create_rotation(qcd_gaugeTransformation* G, qcd_gaugeField* u, qcd_int_4 OrT);
void create_rotationsu2(qcd_gaugeTransformation* Gsu2, qcd_gaugeTransformation* G,\
	qcd_int_4 p, qcd_int_4 q, qcd_int_4 parity, qcd_real_8 overrelax);
void rotate(qcd_gaugeTransformation* Gsu2, qcd_gaugeField* u);
void gaugefixstep(qcd_gaugeField* u, qcd_real_8* current_av, qcd_int_4 OrT, qcd_real_8 overrelax);
void reunitarize_su3(qcd_complex_16 g[3][3]);
void  reunitarize(qcd_gaugeField* u);
void ContinuumField(qcd_gaugeField* A, qcd_gaugeField* u, qcd_int_4 OrT);
qcd_real_8 gaugecheck(qcd_gaugeField* u, qcd_int_4 OrT);
void gaugefix(int gauge_type, qcd_gaugeField* u, qcd_real_8 tol, qcd_int_4 max_iter, qcd_real_8 overrelax);
#endif
