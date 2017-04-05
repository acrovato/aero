#ifndef FPMV1_INTERP_SP_H
#define FPMV1_INTERP_SP_H

#include "Network.h"
#include "Subpanel.h"

Eigen::MatrixXd interp_sp(int idP, Network &bPan, Subpanel &sp,
                   double mu0, double mu1, double mu2, double mu3,
                   double tau0, double tau1, double tau2, double tau3);

#endif //FPMV1_INTERP_SP_H
