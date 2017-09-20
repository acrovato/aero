#ifndef FPMV1_MAP_DERIVATIVES_H
#define FPMV1_MAP_DERIVATIVES_H

#include "Network.h"
#include "Field.h"
#include "Numerical_CST.h"

void map_derivatives(Eigen::MatrixX3d &sGrid, Numerical_CST &numC, Network &bPan, Network &wPan, Field &fPan);

#endif //FPMV1_MAP_DERIVATIVES_H
