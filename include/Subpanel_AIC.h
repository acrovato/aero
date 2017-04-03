//// Subpanel AIC
// Structure containing subpanel AIC matrices (body to field)
// Each vector is of size NSP and contains a matrix of size NS * NF where NSP is the total number of panel to be split
// into sub-panels, NS is the number of sub-panel on a panel and NF is the number of field cell to close to a particular
// panel.

#ifndef FPMV1_SUBPANEL_AIC_H
#define FPMV1_SUBPANEL_AIC_H

#include <vector>
#include <Eigen/Dense>

struct Subpanel_AIC {

    // Cell center
    std::vector <Eigen::MatrixXd> Au;
    std::vector <Eigen::MatrixXd> Av;
    std::vector <Eigen::MatrixXd> Aw;
    std::vector <Eigen::MatrixXd> Bu;
    std::vector <Eigen::MatrixXd> Bv;
    std::vector <Eigen::MatrixXd> Bw;
    // X-line
    std::vector <Eigen::MatrixXd> AuXbwd;
    std::vector <Eigen::MatrixXd> AvXbwd;
    std::vector <Eigen::MatrixXd> AwXbwd;
    std::vector <Eigen::MatrixXd> BuXbwd;
    std::vector <Eigen::MatrixXd> BvXbwd;
    std::vector <Eigen::MatrixXd> BwXbwd;
    std::vector <Eigen::MatrixXd> AuXfwd;
    std::vector <Eigen::MatrixXd> AvXfwd;
    std::vector <Eigen::MatrixXd> AwXfwd;
    std::vector <Eigen::MatrixXd> BuXfwd;
    std::vector <Eigen::MatrixXd> BvXfwd;
    std::vector <Eigen::MatrixXd> BwXfwd;
    // Y-line
    std::vector <Eigen::MatrixXd> AuYbwd;
    std::vector <Eigen::MatrixXd> AvYbwd;
    std::vector <Eigen::MatrixXd> AwYbwd;
    std::vector <Eigen::MatrixXd> BuYbwd;
    std::vector <Eigen::MatrixXd> BvYbwd;
    std::vector <Eigen::MatrixXd> BwYbwd;
    std::vector <Eigen::MatrixXd> AuYfwd;
    std::vector <Eigen::MatrixXd> AvYfwd;
    std::vector <Eigen::MatrixXd> AwYfwd;
    std::vector <Eigen::MatrixXd> BuYfwd;
    std::vector <Eigen::MatrixXd> BvYfwd;
    std::vector <Eigen::MatrixXd> BwYfwd;
    //Z-line
    std::vector <Eigen::MatrixXd> AuZbwd;
    std::vector <Eigen::MatrixXd> AvZbwd;
    std::vector <Eigen::MatrixXd> AwZbwd;
    std::vector <Eigen::MatrixXd> BuZbwd;
    std::vector <Eigen::MatrixXd> BvZbwd;
    std::vector <Eigen::MatrixXd> BwZbwd;
    std::vector <Eigen::MatrixXd> AuZfwd;
    std::vector <Eigen::MatrixXd> AvZfwd;
    std::vector <Eigen::MatrixXd> AwZfwd;
    std::vector <Eigen::MatrixXd> BuZfwd;
    std::vector <Eigen::MatrixXd> BvZfwd;
    std::vector <Eigen::MatrixXd> BwZfwd;
};

#endif //FPMV1_SUBPANEL_AIC_H
