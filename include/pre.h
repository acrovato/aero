#ifndef FPMV1_PRE_H
#define FPMV1_PRE_H

#include "Network.h"
#include "Field.h"

int pre(bool &symY, double &sRef, double &machInf, double &AoA, Eigen::Vector3d &vInf,
        Network &bPan, Network &wPan, Field &fPan);

#endif //FPMV1_PRE_H
