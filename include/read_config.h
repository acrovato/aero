/* Copyright (C) 2018 Adrien Crovato */

#ifndef FPMV1_READ_CONFIG_H
#define FPMV1_READ_CONFIG_H

#include "Numerical_CST.h"
#include "Field.h"
#include "Subpanel.h"

void read_config(std::string path, Numerical_CST &numC, bool &symY, double &sRef, double &machInf, double &AoA,
                 std::array<std::array<double, 3>, 8> &box, Field &fPan, Subpanel &sp);

#endif //FPMV1_READ_CONFIG_H
