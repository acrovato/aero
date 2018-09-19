//// Subpanel AIC
// Structure containing sub-panel AIC matrices (body to field)
// Each vector is of size NSP and contains a matrix of size NS * NF where NSP is the total number of panel to be split
// into sub-panels, NS is the number of sub-panel on a panel and NF is the number of field cell to close to a particular
// panel.

/* Copyright (C) 2018 Adrien Crovato */

#ifndef FPMV1_SUBPANEL_AIC_H
#define FPMV1_SUBPANEL_AIC_H

#include <vector>
#include <Eigen/Dense>

struct Subpanel_AIC {

    // Cell center
    std::vector <Eigen::MatrixXd> A; // doublet matrix
    std::vector <Eigen::MatrixXd> B; // source matrix
};

#endif //FPMV1_SUBPANEL_AIC_H
