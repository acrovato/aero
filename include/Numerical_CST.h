//// Numerical constants structure
// Structure containing numerical constants

/* Copyright (C) 2018 Adrien Crovato */

#ifndef FPMV1_NUMERICAL_CST_H
#define FPMV1_NUMERICAL_CST_H

struct Numerical_CST {
    // Pre-processing
    double TOLB = 1e-3; // geometric tolerance on box (enclosing the geometry) size
    // Solver
    double RRED = 5; // order of magnitude of residual (relative change in field source) reduction
};

#endif //FPMV1_NUMERICAL_CST_H
