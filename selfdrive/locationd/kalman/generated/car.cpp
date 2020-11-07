
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
void err_fun(double *nom_x, double *delta_x, double *out_4234581838440915095) {
   out_4234581838440915095[0] = delta_x[0] + nom_x[0];
   out_4234581838440915095[1] = delta_x[1] + nom_x[1];
   out_4234581838440915095[2] = delta_x[2] + nom_x[2];
   out_4234581838440915095[3] = delta_x[3] + nom_x[3];
   out_4234581838440915095[4] = delta_x[4] + nom_x[4];
   out_4234581838440915095[5] = delta_x[5] + nom_x[5];
   out_4234581838440915095[6] = delta_x[6] + nom_x[6];
   out_4234581838440915095[7] = delta_x[7] + nom_x[7];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_8050171482298030105) {
   out_8050171482298030105[0] = -nom_x[0] + true_x[0];
   out_8050171482298030105[1] = -nom_x[1] + true_x[1];
   out_8050171482298030105[2] = -nom_x[2] + true_x[2];
   out_8050171482298030105[3] = -nom_x[3] + true_x[3];
   out_8050171482298030105[4] = -nom_x[4] + true_x[4];
   out_8050171482298030105[5] = -nom_x[5] + true_x[5];
   out_8050171482298030105[6] = -nom_x[6] + true_x[6];
   out_8050171482298030105[7] = -nom_x[7] + true_x[7];
}
void H_mod_fun(double *state, double *out_4318879525745610550) {
   out_4318879525745610550[0] = 1.0;
   out_4318879525745610550[1] = 0.0;
   out_4318879525745610550[2] = 0.0;
   out_4318879525745610550[3] = 0.0;
   out_4318879525745610550[4] = 0.0;
   out_4318879525745610550[5] = 0.0;
   out_4318879525745610550[6] = 0.0;
   out_4318879525745610550[7] = 0.0;
   out_4318879525745610550[8] = 0.0;
   out_4318879525745610550[9] = 1.0;
   out_4318879525745610550[10] = 0.0;
   out_4318879525745610550[11] = 0.0;
   out_4318879525745610550[12] = 0.0;
   out_4318879525745610550[13] = 0.0;
   out_4318879525745610550[14] = 0.0;
   out_4318879525745610550[15] = 0.0;
   out_4318879525745610550[16] = 0.0;
   out_4318879525745610550[17] = 0.0;
   out_4318879525745610550[18] = 1.0;
   out_4318879525745610550[19] = 0.0;
   out_4318879525745610550[20] = 0.0;
   out_4318879525745610550[21] = 0.0;
   out_4318879525745610550[22] = 0.0;
   out_4318879525745610550[23] = 0.0;
   out_4318879525745610550[24] = 0.0;
   out_4318879525745610550[25] = 0.0;
   out_4318879525745610550[26] = 0.0;
   out_4318879525745610550[27] = 1.0;
   out_4318879525745610550[28] = 0.0;
   out_4318879525745610550[29] = 0.0;
   out_4318879525745610550[30] = 0.0;
   out_4318879525745610550[31] = 0.0;
   out_4318879525745610550[32] = 0.0;
   out_4318879525745610550[33] = 0.0;
   out_4318879525745610550[34] = 0.0;
   out_4318879525745610550[35] = 0.0;
   out_4318879525745610550[36] = 1.0;
   out_4318879525745610550[37] = 0.0;
   out_4318879525745610550[38] = 0.0;
   out_4318879525745610550[39] = 0.0;
   out_4318879525745610550[40] = 0.0;
   out_4318879525745610550[41] = 0.0;
   out_4318879525745610550[42] = 0.0;
   out_4318879525745610550[43] = 0.0;
   out_4318879525745610550[44] = 0.0;
   out_4318879525745610550[45] = 1.0;
   out_4318879525745610550[46] = 0.0;
   out_4318879525745610550[47] = 0.0;
   out_4318879525745610550[48] = 0.0;
   out_4318879525745610550[49] = 0.0;
   out_4318879525745610550[50] = 0.0;
   out_4318879525745610550[51] = 0.0;
   out_4318879525745610550[52] = 0.0;
   out_4318879525745610550[53] = 0.0;
   out_4318879525745610550[54] = 1.0;
   out_4318879525745610550[55] = 0.0;
   out_4318879525745610550[56] = 0.0;
   out_4318879525745610550[57] = 0.0;
   out_4318879525745610550[58] = 0.0;
   out_4318879525745610550[59] = 0.0;
   out_4318879525745610550[60] = 0.0;
   out_4318879525745610550[61] = 0.0;
   out_4318879525745610550[62] = 0.0;
   out_4318879525745610550[63] = 1.0;
}
void f_fun(double *state, double dt, double *out_9213022374902016956) {
   out_9213022374902016956[0] = state[0];
   out_9213022374902016956[1] = state[1];
   out_9213022374902016956[2] = state[2];
   out_9213022374902016956[3] = state[3];
   out_9213022374902016956[4] = state[4];
   out_9213022374902016956[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_9213022374902016956[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_9213022374902016956[7] = state[7];
}
void F_fun(double *state, double dt, double *out_5738822049197680580) {
   out_5738822049197680580[0] = 1;
   out_5738822049197680580[1] = 0;
   out_5738822049197680580[2] = 0;
   out_5738822049197680580[3] = 0;
   out_5738822049197680580[4] = 0;
   out_5738822049197680580[5] = 0;
   out_5738822049197680580[6] = 0;
   out_5738822049197680580[7] = 0;
   out_5738822049197680580[8] = 0;
   out_5738822049197680580[9] = 1;
   out_5738822049197680580[10] = 0;
   out_5738822049197680580[11] = 0;
   out_5738822049197680580[12] = 0;
   out_5738822049197680580[13] = 0;
   out_5738822049197680580[14] = 0;
   out_5738822049197680580[15] = 0;
   out_5738822049197680580[16] = 0;
   out_5738822049197680580[17] = 0;
   out_5738822049197680580[18] = 1;
   out_5738822049197680580[19] = 0;
   out_5738822049197680580[20] = 0;
   out_5738822049197680580[21] = 0;
   out_5738822049197680580[22] = 0;
   out_5738822049197680580[23] = 0;
   out_5738822049197680580[24] = 0;
   out_5738822049197680580[25] = 0;
   out_5738822049197680580[26] = 0;
   out_5738822049197680580[27] = 1;
   out_5738822049197680580[28] = 0;
   out_5738822049197680580[29] = 0;
   out_5738822049197680580[30] = 0;
   out_5738822049197680580[31] = 0;
   out_5738822049197680580[32] = 0;
   out_5738822049197680580[33] = 0;
   out_5738822049197680580[34] = 0;
   out_5738822049197680580[35] = 0;
   out_5738822049197680580[36] = 1;
   out_5738822049197680580[37] = 0;
   out_5738822049197680580[38] = 0;
   out_5738822049197680580[39] = 0;
   out_5738822049197680580[40] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_5738822049197680580[41] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_5738822049197680580[42] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_5738822049197680580[43] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_5738822049197680580[44] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_5738822049197680580[45] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_5738822049197680580[46] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_5738822049197680580[47] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_5738822049197680580[48] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_5738822049197680580[49] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_5738822049197680580[50] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_5738822049197680580[51] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_5738822049197680580[52] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_5738822049197680580[53] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_5738822049197680580[54] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_5738822049197680580[55] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_5738822049197680580[56] = 0;
   out_5738822049197680580[57] = 0;
   out_5738822049197680580[58] = 0;
   out_5738822049197680580[59] = 0;
   out_5738822049197680580[60] = 0;
   out_5738822049197680580[61] = 0;
   out_5738822049197680580[62] = 0;
   out_5738822049197680580[63] = 1;
}
void h_25(double *state, double *unused, double *out_3542551237994567859) {
   out_3542551237994567859[0] = state[6];
}
void H_25(double *state, double *unused, double *out_4859059977381893746) {
   out_4859059977381893746[0] = 0;
   out_4859059977381893746[1] = 0;
   out_4859059977381893746[2] = 0;
   out_4859059977381893746[3] = 0;
   out_4859059977381893746[4] = 0;
   out_4859059977381893746[5] = 0;
   out_4859059977381893746[6] = 1;
   out_4859059977381893746[7] = 0;
}
void h_24(double *state, double *unused, double *out_6760447318376270888) {
   out_6760447318376270888[0] = state[4];
   out_6760447318376270888[1] = state[5];
}
void H_24(double *state, double *unused, double *out_6255321910726500370) {
   out_6255321910726500370[0] = 0;
   out_6255321910726500370[1] = 0;
   out_6255321910726500370[2] = 0;
   out_6255321910726500370[3] = 0;
   out_6255321910726500370[4] = 1;
   out_6255321910726500370[5] = 0;
   out_6255321910726500370[6] = 0;
   out_6255321910726500370[7] = 0;
   out_6255321910726500370[8] = 0;
   out_6255321910726500370[9] = 0;
   out_6255321910726500370[10] = 0;
   out_6255321910726500370[11] = 0;
   out_6255321910726500370[12] = 0;
   out_6255321910726500370[13] = 1;
   out_6255321910726500370[14] = 0;
   out_6255321910726500370[15] = 0;
}
void h_30(double *state, double *unused, double *out_5489264300762674993) {
   out_5489264300762674993[0] = state[4];
}
void H_30(double *state, double *unused, double *out_3675621562500390728) {
   out_3675621562500390728[0] = 0;
   out_3675621562500390728[1] = 0;
   out_3675621562500390728[2] = 0;
   out_3675621562500390728[3] = 0;
   out_3675621562500390728[4] = 1;
   out_3675621562500390728[5] = 0;
   out_3675621562500390728[6] = 0;
   out_3675621562500390728[7] = 0;
}
void h_26(double *state, double *unused, double *out_2909342249642507268) {
   out_2909342249642507268[0] = state[7];
}
void H_26(double *state, double *unused, double *out_5281939879324821562) {
   out_5281939879324821562[0] = 0;
   out_5281939879324821562[1] = 0;
   out_5281939879324821562[2] = 0;
   out_5281939879324821562[3] = 0;
   out_5281939879324821562[4] = 0;
   out_5281939879324821562[5] = 0;
   out_5281939879324821562[6] = 0;
   out_5281939879324821562[7] = 1;
}
void h_27(double *state, double *unused, double *out_7289276073878845840) {
   out_7289276073878845840[0] = state[3];
}
void H_27(double *state, double *unused, double *out_927550513521344618) {
   out_927550513521344618[0] = 0;
   out_927550513521344618[1] = 0;
   out_927550513521344618[2] = 0;
   out_927550513521344618[3] = 1;
   out_927550513521344618[4] = 0;
   out_927550513521344618[5] = 0;
   out_927550513521344618[6] = 0;
   out_927550513521344618[7] = 0;
}
void h_29(double *state, double *unused, double *out_8354707923849902866) {
   out_8354707923849902866[0] = state[1];
}
void H_29(double *state, double *unused, double *out_9122555981703708322) {
   out_9122555981703708322[0] = 0;
   out_9122555981703708322[1] = 1;
   out_9122555981703708322[2] = 0;
   out_9122555981703708322[3] = 0;
   out_9122555981703708322[4] = 0;
   out_9122555981703708322[5] = 0;
   out_9122555981703708322[6] = 0;
   out_9122555981703708322[7] = 0;
}
void h_28(double *state, double *unused, double *out_4083533142336211128) {
   out_4083533142336211128[0] = state[5];
   out_4083533142336211128[1] = state[6];
}
void H_28(double *state, double *unused, double *out_6056233691103076300) {
   out_6056233691103076300[0] = 0;
   out_6056233691103076300[1] = 0;
   out_6056233691103076300[2] = 0;
   out_6056233691103076300[3] = 0;
   out_6056233691103076300[4] = 0;
   out_6056233691103076300[5] = 1;
   out_6056233691103076300[6] = 0;
   out_6056233691103076300[7] = 0;
   out_6056233691103076300[8] = 0;
   out_6056233691103076300[9] = 0;
   out_6056233691103076300[10] = 0;
   out_6056233691103076300[11] = 0;
   out_6056233691103076300[12] = 0;
   out_6056233691103076300[13] = 0;
   out_6056233691103076300[14] = 1;
   out_6056233691103076300[15] = 0;
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
