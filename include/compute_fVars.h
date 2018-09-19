/* Copyright (C) 2018 Adrien Crovato */

#ifndef FPMV1_COMPUTE_FVARS_H
#define FPMV1_COMPUTE_FVARS_H

#include "Network.h"
#include "Field.h"
#include "Subpanel.h"
#include "Body_AIC.h"
#include "Field2field_AIC.h"
#include "Subpanel_AIC.h"

void compute_fVars(double Minf, Eigen::Vector3d &vInf, Network &bPan, Field &fPan, Subpanel &sp,
                       Body_AIC &b2fAIC, Field2field_AIC &f2fAIC, Subpanel_AIC &spAIC);

#endif //FPMV1_COMPUTE_FVARS_H
