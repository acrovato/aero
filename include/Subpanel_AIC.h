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
};

#endif //FPMV1_SUBPANEL_AIC_H
