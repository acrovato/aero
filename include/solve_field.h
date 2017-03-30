#ifndef FPMV1_COMP_FIELD_H
#define FPMV1_COMP_FIELD_H

#include "Network.h"
#include "Field.h"
#include "Minigrid.h"
#include "Field_AIC.h"
#include "Minigrid_AIC.h"
#include "Body2field_AIC.h"

void solve_field(double Minf, Eigen::Vector3d &vInf, Network &bPan, Field &fPan, Minigrid &mgVar, Eigen::MatrixX3d &dRho,
                 Body2field_AIC &b2fAIC, Field_AIC &f2fAIC, Minigrid_AIC &mgAIC);

#endif //FPMV1_COMP_FIELD_H
