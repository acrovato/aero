//// Body to body/field AIC structure
// Structure containing body to body/field AIC matrices

/* Copyright (C) 2018 Adrien Crovato */

#ifndef FPMV1_BODY_AIC_H
#define FPMV1_BODY_AIC_H

#include <Eigen/Dense>

struct Body_AIC {
    Eigen::MatrixXd A; // doublet AIC matrix
    Eigen::MatrixXd B; // source AIC matrix
};

#endif //FPMV1_BODY_AIC_H
