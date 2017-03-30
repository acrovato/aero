#ifndef FPMV1_MINIGRID_AIC_H
#define FPMV1_MINIGRID_AIC_H

#include <Eigen/Dense>

struct Minigrid_AIC {
    // X-line
    Eigen::MatrixXd AuXbwd;
    Eigen::MatrixXd AvXbwd;
    Eigen::MatrixXd AwXbwd;
    Eigen::MatrixXd BuXbwd;
    Eigen::MatrixXd BvXbwd;
    Eigen::MatrixXd BwXbwd;
    Eigen::MatrixXd AuXfwd;
    Eigen::MatrixXd AvXfwd;
    Eigen::MatrixXd AwXfwd;
    Eigen::MatrixXd BuXfwd;
    Eigen::MatrixXd BvXfwd;
    Eigen::MatrixXd BwXfwd;
    Eigen::MatrixXd CuXbwd;
    Eigen::MatrixXd CvXbwd;
    Eigen::MatrixXd CwXbwd;
    Eigen::MatrixXd CuXfwd;
    Eigen::MatrixXd CvXfwd;
    Eigen::MatrixXd CwXfwd;
    // Y-line
    Eigen::MatrixXd AuYbwd;
    Eigen::MatrixXd AvYbwd;
    Eigen::MatrixXd AwYbwd;
    Eigen::MatrixXd BuYbwd;
    Eigen::MatrixXd BvYbwd;
    Eigen::MatrixXd BwYbwd;
    Eigen::MatrixXd AuYfwd;
    Eigen::MatrixXd AvYfwd;
    Eigen::MatrixXd AwYfwd;
    Eigen::MatrixXd BuYfwd;
    Eigen::MatrixXd BvYfwd;
    Eigen::MatrixXd BwYfwd;
    Eigen::MatrixXd CuYbwd;
    Eigen::MatrixXd CvYbwd;
    Eigen::MatrixXd CwYbwd;
    Eigen::MatrixXd CuYfwd;
    Eigen::MatrixXd CvYfwd;
    Eigen::MatrixXd CwYfwd;
    // Z-line
    Eigen::MatrixXd AuZbwd;
    Eigen::MatrixXd AvZbwd;
    Eigen::MatrixXd AwZbwd;
    Eigen::MatrixXd BuZbwd;
    Eigen::MatrixXd BvZbwd;
    Eigen::MatrixXd BwZbwd;
    Eigen::MatrixXd AuZfwd;
    Eigen::MatrixXd AvZfwd;
    Eigen::MatrixXd AwZfwd;
    Eigen::MatrixXd BuZfwd;
    Eigen::MatrixXd BvZfwd;
    Eigen::MatrixXd BwZfwd;
    Eigen::MatrixXd CuZbwd;
    Eigen::MatrixXd CvZbwd;
    Eigen::MatrixXd CwZbwd;
    Eigen::MatrixXd CuZfwd;
    Eigen::MatrixXd CvZfwd;
    Eigen::MatrixXd CwZfwd;
};

#endif //FPMV1_MINIGRID_AIC_H

