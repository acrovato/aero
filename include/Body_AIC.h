#ifndef FPMV1_BODY_AIC_H
#define FPMV1_BODY_AIC_H

#include <Eigen/Dense>

struct Body_AIC {
    Eigen::MatrixXd A;
    Eigen::MatrixXd B;
};

#endif //FPMV1_BODY_AIC_H
