#ifndef FPMV1_READ_CONFIG_H
#define FPMV1_READ_CONFIG_H

#include "Field.h"

void read_config(std::string path, bool &symY, double &sRef, double &machInf, double &AoA,
                 std::array<std::array<double, 3>, 8> &box, Field &fPan);

#endif //FPMV1_READ_CONFIG_H
