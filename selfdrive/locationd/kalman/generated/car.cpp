
extern "C"{

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}

}
extern "C" {
#include <math.h>
/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_6109257089592663609) {
   out_6109257089592663609[0] = delta_x[0] + nom_x[0];
   out_6109257089592663609[1] = delta_x[1] + nom_x[1];
   out_6109257089592663609[2] = delta_x[2] + nom_x[2];
   out_6109257089592663609[3] = delta_x[3] + nom_x[3];
   out_6109257089592663609[4] = delta_x[4] + nom_x[4];
   out_6109257089592663609[5] = delta_x[5] + nom_x[5];
   out_6109257089592663609[6] = delta_x[6] + nom_x[6];
   out_6109257089592663609[7] = delta_x[7] + nom_x[7];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_2675833930546276745) {
   out_2675833930546276745[0] = -nom_x[0] + true_x[0];
   out_2675833930546276745[1] = -nom_x[1] + true_x[1];
   out_2675833930546276745[2] = -nom_x[2] + true_x[2];
   out_2675833930546276745[3] = -nom_x[3] + true_x[3];
   out_2675833930546276745[4] = -nom_x[4] + true_x[4];
   out_2675833930546276745[5] = -nom_x[5] + true_x[5];
   out_2675833930546276745[6] = -nom_x[6] + true_x[6];
   out_2675833930546276745[7] = -nom_x[7] + true_x[7];
}
void H_mod_fun(double *state, double *out_7193349714261281782) {
   out_7193349714261281782[0] = 1.0;
   out_7193349714261281782[1] = 0.0;
   out_7193349714261281782[2] = 0.0;
   out_7193349714261281782[3] = 0.0;
   out_7193349714261281782[4] = 0.0;
   out_7193349714261281782[5] = 0.0;
   out_7193349714261281782[6] = 0.0;
   out_7193349714261281782[7] = 0.0;
   out_7193349714261281782[8] = 0.0;
   out_7193349714261281782[9] = 1.0;
   out_7193349714261281782[10] = 0.0;
   out_7193349714261281782[11] = 0.0;
   out_7193349714261281782[12] = 0.0;
   out_7193349714261281782[13] = 0.0;
   out_7193349714261281782[14] = 0.0;
   out_7193349714261281782[15] = 0.0;
   out_7193349714261281782[16] = 0.0;
   out_7193349714261281782[17] = 0.0;
   out_7193349714261281782[18] = 1.0;
   out_7193349714261281782[19] = 0.0;
   out_7193349714261281782[20] = 0.0;
   out_7193349714261281782[21] = 0.0;
   out_7193349714261281782[22] = 0.0;
   out_7193349714261281782[23] = 0.0;
   out_7193349714261281782[24] = 0.0;
   out_7193349714261281782[25] = 0.0;
   out_7193349714261281782[26] = 0.0;
   out_7193349714261281782[27] = 1.0;
   out_7193349714261281782[28] = 0.0;
   out_7193349714261281782[29] = 0.0;
   out_7193349714261281782[30] = 0.0;
   out_7193349714261281782[31] = 0.0;
   out_7193349714261281782[32] = 0.0;
   out_7193349714261281782[33] = 0.0;
   out_7193349714261281782[34] = 0.0;
   out_7193349714261281782[35] = 0.0;
   out_7193349714261281782[36] = 1.0;
   out_7193349714261281782[37] = 0.0;
   out_7193349714261281782[38] = 0.0;
   out_7193349714261281782[39] = 0.0;
   out_7193349714261281782[40] = 0.0;
   out_7193349714261281782[41] = 0.0;
   out_7193349714261281782[42] = 0.0;
   out_7193349714261281782[43] = 0.0;
   out_7193349714261281782[44] = 0.0;
   out_7193349714261281782[45] = 1.0;
   out_7193349714261281782[46] = 0.0;
   out_7193349714261281782[47] = 0.0;
   out_7193349714261281782[48] = 0.0;
   out_7193349714261281782[49] = 0.0;
   out_7193349714261281782[50] = 0.0;
   out_7193349714261281782[51] = 0.0;
   out_7193349714261281782[52] = 0.0;
   out_7193349714261281782[53] = 0.0;
   out_7193349714261281782[54] = 1.0;
   out_7193349714261281782[55] = 0.0;
   out_7193349714261281782[56] = 0.0;
   out_7193349714261281782[57] = 0.0;
   out_7193349714261281782[58] = 0.0;
   out_7193349714261281782[59] = 0.0;
   out_7193349714261281782[60] = 0.0;
   out_7193349714261281782[61] = 0.0;
   out_7193349714261281782[62] = 0.0;
   out_7193349714261281782[63] = 1.0;
}
void f_fun(double *state, double dt, double *out_4864488243241977304) {
   out_4864488243241977304[0] = state[0];
   out_4864488243241977304[1] = state[1];
   out_4864488243241977304[2] = state[2];
   out_4864488243241977304[3] = state[3];
   out_4864488243241977304[4] = state[4];
   out_4864488243241977304[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_4864488243241977304[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_4864488243241977304[7] = state[7];
}
void F_fun(double *state, double dt, double *out_2328359433604071268) {
   out_2328359433604071268[0] = 1;
   out_2328359433604071268[1] = 0;
   out_2328359433604071268[2] = 0;
   out_2328359433604071268[3] = 0;
   out_2328359433604071268[4] = 0;
   out_2328359433604071268[5] = 0;
   out_2328359433604071268[6] = 0;
   out_2328359433604071268[7] = 0;
   out_2328359433604071268[8] = 0;
   out_2328359433604071268[9] = 1;
   out_2328359433604071268[10] = 0;
   out_2328359433604071268[11] = 0;
   out_2328359433604071268[12] = 0;
   out_2328359433604071268[13] = 0;
   out_2328359433604071268[14] = 0;
   out_2328359433604071268[15] = 0;
   out_2328359433604071268[16] = 0;
   out_2328359433604071268[17] = 0;
   out_2328359433604071268[18] = 1;
   out_2328359433604071268[19] = 0;
   out_2328359433604071268[20] = 0;
   out_2328359433604071268[21] = 0;
   out_2328359433604071268[22] = 0;
   out_2328359433604071268[23] = 0;
   out_2328359433604071268[24] = 0;
   out_2328359433604071268[25] = 0;
   out_2328359433604071268[26] = 0;
   out_2328359433604071268[27] = 1;
   out_2328359433604071268[28] = 0;
   out_2328359433604071268[29] = 0;
   out_2328359433604071268[30] = 0;
   out_2328359433604071268[31] = 0;
   out_2328359433604071268[32] = 0;
   out_2328359433604071268[33] = 0;
   out_2328359433604071268[34] = 0;
   out_2328359433604071268[35] = 0;
   out_2328359433604071268[36] = 1;
   out_2328359433604071268[37] = 0;
   out_2328359433604071268[38] = 0;
   out_2328359433604071268[39] = 0;
   out_2328359433604071268[40] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_2328359433604071268[41] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_2328359433604071268[42] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_2328359433604071268[43] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_2328359433604071268[44] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_2328359433604071268[45] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_2328359433604071268[46] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_2328359433604071268[47] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_2328359433604071268[48] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_2328359433604071268[49] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_2328359433604071268[50] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_2328359433604071268[51] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_2328359433604071268[52] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_2328359433604071268[53] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_2328359433604071268[54] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_2328359433604071268[55] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_2328359433604071268[56] = 0;
   out_2328359433604071268[57] = 0;
   out_2328359433604071268[58] = 0;
   out_2328359433604071268[59] = 0;
   out_2328359433604071268[60] = 0;
   out_2328359433604071268[61] = 0;
   out_2328359433604071268[62] = 0;
   out_2328359433604071268[63] = 1;
}
void h_25(double *state, double *unused, double *out_312706834215263469) {
   out_312706834215263469[0] = state[6];
}
void H_25(double *state, double *unused, double *out_7812137926698848974) {
   out_7812137926698848974[0] = 0;
   out_7812137926698848974[1] = 0;
   out_7812137926698848974[2] = 0;
   out_7812137926698848974[3] = 0;
   out_7812137926698848974[4] = 0;
   out_7812137926698848974[5] = 0;
   out_7812137926698848974[6] = 1;
   out_7812137926698848974[7] = 0;
}
void h_24(double *state, double *unused, double *out_9026545562781906584) {
   out_9026545562781906584[0] = state[4];
   out_9026545562781906584[1] = state[5];
}
void H_24(double *state, double *unused, double *out_5342461122073738670) {
   out_5342461122073738670[0] = 0;
   out_5342461122073738670[1] = 0;
   out_5342461122073738670[2] = 0;
   out_5342461122073738670[3] = 0;
   out_5342461122073738670[4] = 1;
   out_5342461122073738670[5] = 0;
   out_5342461122073738670[6] = 0;
   out_5342461122073738670[7] = 0;
   out_5342461122073738670[8] = 0;
   out_5342461122073738670[9] = 0;
   out_5342461122073738670[10] = 0;
   out_5342461122073738670[11] = 0;
   out_5342461122073738670[12] = 0;
   out_5342461122073738670[13] = 1;
   out_5342461122073738670[14] = 0;
   out_5342461122073738670[15] = 0;
}
void h_30(double *state, double *unused, double *out_7803559022344293009) {
   out_7803559022344293009[0] = state[4];
}
void H_30(double *state, double *unused, double *out_7655460888536193528) {
   out_7655460888536193528[0] = 0;
   out_7655460888536193528[1] = 0;
   out_7655460888536193528[2] = 0;
   out_7655460888536193528[3] = 0;
   out_7655460888536193528[4] = 1;
   out_7655460888536193528[5] = 0;
   out_7655460888536193528[6] = 0;
   out_7655460888536193528[7] = 0;
}
void h_26(double *state, double *unused, double *out_8657089636458053340) {
   out_8657089636458053340[0] = state[7];
}
void H_26(double *state, double *unused, double *out_3624683376761947642) {
   out_3624683376761947642[0] = 0;
   out_3624683376761947642[1] = 0;
   out_3624683376761947642[2] = 0;
   out_3624683376761947642[3] = 0;
   out_3624683376761947642[4] = 0;
   out_3624683376761947642[5] = 0;
   out_3624683376761947642[6] = 0;
   out_3624683376761947642[7] = 1;
}
void h_27(double *state, double *unused, double *out_8369073581914504400) {
   out_8369073581914504400[0] = state[3];
}
void H_27(double *state, double *unused, double *out_4459501234083275434) {
   out_4459501234083275434[0] = 0;
   out_4459501234083275434[1] = 0;
   out_4459501234083275434[2] = 0;
   out_4459501234083275434[3] = 1;
   out_4459501234083275434[4] = 0;
   out_4459501234083275434[5] = 0;
   out_4459501234083275434[6] = 0;
   out_4459501234083275434[7] = 0;
}
void h_29(double *state, double *unused, double *out_6224705531520762958) {
   out_6224705531520762958[0] = state[1];
}
void H_29(double *state, double *unused, double *out_3575685844487292958) {
   out_3575685844487292958[0] = 0;
   out_3575685844487292958[1] = 1;
   out_3575685844487292958[2] = 0;
   out_3575685844487292958[3] = 0;
   out_3575685844487292958[4] = 0;
   out_3575685844487292958[5] = 0;
   out_3575685844487292958[6] = 0;
   out_3575685844487292958[7] = 0;
}
void h_28(double *state, double *unused, double *out_1705069466177101320) {
   out_1705069466177101320[0] = state[5];
   out_1705069466177101320[1] = state[6];
}
void H_28(double *state, double *unused, double *out_1289310505460256244) {
   out_1289310505460256244[0] = 0;
   out_1289310505460256244[1] = 0;
   out_1289310505460256244[2] = 0;
   out_1289310505460256244[3] = 0;
   out_1289310505460256244[4] = 0;
   out_1289310505460256244[5] = 1;
   out_1289310505460256244[6] = 0;
   out_1289310505460256244[7] = 0;
   out_1289310505460256244[8] = 0;
   out_1289310505460256244[9] = 0;
   out_1289310505460256244[10] = 0;
   out_1289310505460256244[11] = 0;
   out_1289310505460256244[12] = 0;
   out_1289310505460256244[13] = 0;
   out_1289310505460256244[14] = 1;
   out_1289310505460256244[15] = 0;
}
}

extern "C"{
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
}

#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;
  
  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);
  
  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H); 
  
  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();
   

    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;
  
  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);
 
  // update cov 
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}



extern "C"{

      void update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
      }
    
      void update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
      }
    
      void update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
      }
    
      void update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
      }
    
      void update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
      }
    
      void update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
      }
    
      void update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
      }
    
}
