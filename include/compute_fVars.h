#ifndef FPMV1_COMPUTE_FVARS_H
#define FPMV1_COMPUTE_FVARS_H

#include "Network.h"
#include "Field.h"
#include "Minigrid.h"
#include "Body2field_AIC.h"
#include "Field_AIC.h"
#include "Minigrid_AIC.h"

void compute_fVars(double Minf, Eigen::Vector3d &vInf, Network &bPan, Field &fPan, Minigrid &mgVar, Eigen::MatrixX3d &dRho,
                       Body2field_AIC &b2fAIC, Field_AIC &f2fAIC, Minigrid_AIC &mgAIC);

#endif //FPMV1_COMPUTE_FVARS_H
