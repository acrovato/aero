//// Field to body AIC structure
// Structure containing field to body AIC matrices

/* Copyright (C) 2018 Adrien Crovato */

#ifndef FPMV1_FIELD2BODY_AIC_H
#define FPMV1_FIELD2BODY_AIC_H

#include <Eigen/Dense>

struct Field2body_AIC {
    Eigen::MatrixXd Cu; // source matrix for x-component of velocity
    Eigen::MatrixXd Cv; // source matrix for y-component of velocity
    Eigen::MatrixXd Cw; // source matrix for z-component of velocity
};

#endif //FPMV1_FIELD2BODY_AIC_H
