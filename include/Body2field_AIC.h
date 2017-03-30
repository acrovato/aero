#ifndef FPMV1_BODY2FIELD_AIC_H
#define FPMV1_BODY2FIELD_AIC_H

#include <Eigen/Dense>

struct Body2field_AIC {
    Eigen::MatrixXd Au;
    Eigen::MatrixXd Av;
    Eigen::MatrixXd Aw;
    Eigen::MatrixXd Bu;
    Eigen::MatrixXd Bv;
    Eigen::MatrixXd Bw;
};

#endif //FPMV1_BODY2FIELD_AIC_H
