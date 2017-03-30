#ifndef FPMV1_SOLVER_H
#define FPMV1_SOLVER_H

#include "Network.h"
#include "Field.h"

int solver(bool symY, double sRef, double alpha, Eigen::Vector3d &vInf, double Minf,
           Network &bPan, Network &wPan, Field &fPan, double &cL, double &cD);

#endif //FPMV1_SOLVER_H
