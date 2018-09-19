/* Copyright (C) 2018 Adrien Crovato */

#ifndef FPMV1_CREATE_FIELD_H
#define FPMV1_CREATE_FIELD_H

#include "Field.h"
#include "Numerical_CST.h"

void create_field(std::array<std::array<double, 3>, 8> &box, Field &fPan);

#endif //FPMV1_CREATE_FIELD_H
