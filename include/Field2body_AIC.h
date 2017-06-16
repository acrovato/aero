#ifndef FPMV1_FIELD2BODY_AIC_H
#define FPMV1_FIELD2BODY_AIC_H

#include <Eigen/Dense>

struct Field2body_AIC {
    Eigen::MatrixXd Cu;
    Eigen::MatrixXd Cv;
    Eigen::MatrixXd Cw;
};

#endif //FPMV1_FIELD2BODY_AIC_H
