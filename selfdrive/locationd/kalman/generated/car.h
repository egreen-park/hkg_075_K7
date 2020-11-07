/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_4234581838440915095);
void inv_err_fun(double *nom_x, double *true_x, double *out_8050171482298030105);
void H_mod_fun(double *state, double *out_4318879525745610550);
void f_fun(double *state, double dt, double *out_9213022374902016956);
void F_fun(double *state, double dt, double *out_5738822049197680580);
void h_25(double *state, double *unused, double *out_3542551237994567859);
void H_25(double *state, double *unused, double *out_4859059977381893746);
void h_24(double *state, double *unused, double *out_6760447318376270888);
void H_24(double *state, double *unused, double *out_6255321910726500370);
void h_30(double *state, double *unused, double *out_5489264300762674993);
void H_30(double *state, double *unused, double *out_3675621562500390728);
void h_26(double *state, double *unused, double *out_2909342249642507268);
void H_26(double *state, double *unused, double *out_5281939879324821562);
void h_27(double *state, double *unused, double *out_7289276073878845840);
void H_27(double *state, double *unused, double *out_927550513521344618);
void h_29(double *state, double *unused, double *out_8354707923849902866);
void H_29(double *state, double *unused, double *out_9122555981703708322);
void h_28(double *state, double *unused, double *out_4083533142336211128);
void H_28(double *state, double *unused, double *out_6056233691103076300);
#define DIM 8
#define EDIM 8
#define MEDIM 8
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_25 = 3.841459;
void update_25(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_24 = 5.991465;
void update_24(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_30 = 3.841459;
void update_30(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_26 = 3.841459;
void update_26(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_27 = 3.841459;
void update_27(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_29 = 3.841459;
void update_29(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_28 = 5.991465;
void update_28(double *, double *, double *, double *, double *);
void set_mass(double x);

void set_rotational_inertia(double x);

void set_center_to_front(double x);

void set_center_to_rear(double x);

void set_stiffness_front(double x);

void set_stiffness_rear(double x);
