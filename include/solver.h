#ifndef FPMV1_SOLVER_H
#define FPMV1_SOLVER_H

#include "Network.h"
#include "Field.h"
#include "Numerical_CST.h"
#include "Subpanel.h"

int solver(Numerical_CST &numC, bool symY, double sRef, double alpha, Eigen::Vector3d &vInf, double Minf,
           Network &bPan, Network &wPan, Field &fPan, Subpanel &sp, double &cL, double &cD);

#endif //FPMV1_SOLVER_H
