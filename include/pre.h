#ifndef FPMV1_PRE_H
#define FPMV1_PRE_H

#include "Network.h"
#include "Field.h"
#include "Numerical_CST.h"
#include "Subpanel.h"

int pre(char *argv[], Numerical_CST &numC, bool &symY, double &sRef, double &machInf, double &AoA, Eigen::Vector3d &vInf,
        Network &bPan, Network &wPan, Field &fPan, Subpanel &sp);

#endif //FPMV1_PRE_H
