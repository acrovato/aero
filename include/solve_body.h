/* Copyright (C) 2018 Adrien Crovato */

#ifndef FPMV1_COMP_BODY_H
#define FPMV1_COMP_BODY_H

#include "Network.h"
#include "Body_AIC.h"

void solve_body(Eigen::Vector3d &vInf, Eigen::VectorXd &RHS, Eigen::MatrixX3d &vSigma,
                Network &bPan, Body_AIC &b2bAIC);

#endif //FPMV1_COMP_BODY_H
