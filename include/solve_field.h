#ifndef FPMV1_COMP_FIELD_H
#define FPMV1_COMP_FIELD_H

#include "Network.h"
#include "Field.h"
#include "Subpanel.h"
#include "Field_AIC.h"
#include "Body2field_AIC.h"
#include "Subpanel_AIC.h"

void solve_field(double Minf, Eigen::Vector3d &vInf, Network &bPan, Field &fPan, Subpanel &sp,
                 Body2field_AIC &b2fAIC, Field_AIC &f2fAIC, Subpanel_AIC &spAIC);

#endif //FPMV1_COMP_FIELD_H
