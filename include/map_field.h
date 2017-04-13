#ifndef FPMV1_MAP_FIELD_H
#define FPMV1_MAP_FIELD_H

#include "Field.h"
#include "Network.h"
#include "Numerical_CST.h"

void map_field(Eigen::MatrixX3d &sGrid, Numerical_CST &numC, Network &bPan, Field &fPan);

#endif //FPMV1_MAP_FIELD_H
