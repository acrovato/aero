//// Numerical constants structure
// Structure containing numerical constants

#ifndef FPMV1_NUMERICAL_CST_H
#define FPMV1_NUMERICAL_CST_H

struct Numerical_CST {
    // Pre-processing
    double TOLB = 1e-3; // geometric tolerance on box (enclosing the geometry) size
    double TOLS = 1e-6; // geometric tolerance on distance between field singularity (field panel center) and surface
    // Solver
    double RRED = 5; // order of magnitude of residual (relative change in field source) reduction
};

#endif //FPMV1_NUMERICAL_CST_H
