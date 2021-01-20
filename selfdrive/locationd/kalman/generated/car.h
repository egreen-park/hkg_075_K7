/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_6109257089592663609);
void inv_err_fun(double *nom_x, double *true_x, double *out_2675833930546276745);
void H_mod_fun(double *state, double *out_7193349714261281782);
void f_fun(double *state, double dt, double *out_4864488243241977304);
void F_fun(double *state, double dt, double *out_2328359433604071268);
void h_25(double *state, double *unused, double *out_312706834215263469);
void H_25(double *state, double *unused, double *out_7812137926698848974);
void h_24(double *state, double *unused, double *out_9026545562781906584);
void H_24(double *state, double *unused, double *out_5342461122073738670);
void h_30(double *state, double *unused, double *out_7803559022344293009);
void H_30(double *state, double *unused, double *out_7655460888536193528);
void h_26(double *state, double *unused, double *out_8657089636458053340);
void H_26(double *state, double *unused, double *out_3624683376761947642);
void h_27(double *state, double *unused, double *out_8369073581914504400);
void H_27(double *state, double *unused, double *out_4459501234083275434);
void h_29(double *state, double *unused, double *out_6224705531520762958);
void H_29(double *state, double *unused, double *out_3575685844487292958);
void h_28(double *state, double *unused, double *out_1705069466177101320);
void H_28(double *state, double *unused, double *out_1289310505460256244);
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
