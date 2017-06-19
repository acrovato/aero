//// Field to field AIC structure
// Structure containing field to field AIC matrices

#ifndef FPMV1_FIELD2FIELD_AIC_H
#define FPMV1_FIELD2FIELD_AIC_H

#include <Eigen/Dense>

struct Field2field_AIC {
    Eigen::MatrixXd C; // source matrix
};

#endif //FPMV1_FIELD2FIELD_AIC_H
