/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_8641944843224565383);
void inv_err_fun(double *nom_x, double *true_x, double *out_3651955592106735463);
void H_mod_fun(double *state, double *out_208794138962305086);
void f_fun(double *state, double dt, double *out_8512489912926436601);
void F_fun(double *state, double dt, double *out_146540316369868626);
void h_3(double *state, double *unused, double *out_6003100974745567033);
void H_3(double *state, double *unused, double *out_818615590300771559);
void h_4(double *state, double *unused, double *out_2539820227968139182);
void H_4(double *state, double *unused, double *out_3328432233438209070);
void h_9(double *state, double *unused, double *out_6195079780404799197);
void H_9(double *state, double *unused, double *out_7559482845038024208);
void h_10(double *state, double *unused, double *out_8265722886709985110);
void H_10(double *state, double *unused, double *out_2196057604427071322);
void h_12(double *state, double *unused, double *out_164557912020723271);
void H_12(double *state, double *unused, double *out_2275624875751028640);
void h_13(double *state, double *unused, double *out_8228533971597432968);
void H_13(double *state, double *unused, double *out_5594817917920712136);
void h_14(double *state, double *unused, double *out_6195079780404799197);
void H_14(double *state, double *unused, double *out_7559482845038024208);
void h_19(double *state, double *unused, double *out_6230802499590243163);
void H_19(double *state, double *unused, double *out_3433979606982675220);
#define DIM 23
#define EDIM 22
#define MEDIM 22
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_3 = 3.841459;
void update_3(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_4 = 7.814728;
void update_4(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_9 = 7.814728;
void update_9(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_10 = 7.814728;
void update_10(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_12 = 7.814728;
void update_12(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_13 = 7.814728;
void update_13(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_14 = 7.814728;
void update_14(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_19 = 7.814728;
void update_19(double *, double *, double *, double *, double *);