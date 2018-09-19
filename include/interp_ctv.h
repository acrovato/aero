/* Copyright (C) 2018 Adrien Crovato */

#ifndef FPMV1_INTERP_CTV_H
#define FPMV1_INTERP_CTV_H

#include "Network.h"

Eigen::MatrixXd interp_ctv(int idG, int idC, int idS, Network &bPan);

#endif //FPMV1_INTERP_CTV_H
