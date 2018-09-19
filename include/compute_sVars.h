/* Copyright (C) 2018 Adrien Crovato */

#ifndef FPMV1_COMPUTE_SVARS_H
#define FPMV1_COMPUTE_SVARS_H

#include "Network.h"

void compute_sVars(bool symY, double sRef, double alpha, double Minf, Eigen::Vector3d &vInf,
                Eigen::MatrixX3d &vSigma, Network &bPan, double &cL, double &cD);

#endif //FPMV1_COMPUTE_SQ_H
