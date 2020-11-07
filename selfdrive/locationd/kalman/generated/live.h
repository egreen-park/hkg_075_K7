/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_98376307135833035);
void inv_err_fun(double *nom_x, double *true_x, double *out_8751769822747665453);
void H_mod_fun(double *state, double *out_5518803763119843260);
void f_fun(double *state, double dt, double *out_2683054986129145345);
void F_fun(double *state, double dt, double *out_8068501029736161626);
void h_3(double *state, double *unused, double *out_1115124464371072775);
void H_3(double *state, double *unused, double *out_8521270738708131793);
void h_4(double *state, double *unused, double *out_3047349443930758882);
void H_4(double *state, double *unused, double *out_4795726417735365690);
void h_9(double *state, double *unused, double *out_5708378332327032327);
void H_9(double *state, double *unused, double *out_1737754398311667290);
void h_10(double *state, double *unused, double *out_1980395709464591812);
void H_10(double *state, double *unused, double *out_5504191171440146950);
void h_12(double *state, double *unused, double *out_4129799582437393595);
void H_12(double *state, double *unused, double *out_4308003494150801738);
void h_13(double *state, double *unused, double *out_7336077623459016708);
void H_13(double *state, double *unused, double *out_3105401458294534158);
void h_14(double *state, double *unused, double *out_5708378332327032327);
void H_14(double *state, double *unused, double *out_1737754398311667290);
void h_19(double *state, double *unused, double *out_2107634465980273759);
void H_19(double *state, double *unused, double *out_7523987201661649650);
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