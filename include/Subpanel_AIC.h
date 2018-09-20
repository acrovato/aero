/* 
// Copyright 2018 University of Liege
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// 
//     http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
// Authors:
// - Adrien Crovato
*/

//// Subpanel AIC
// Structure containing sub-panel AIC matrices (body to field)
// Each vector is of size NSP and contains a matrix of size NS * NF where NSP is the total number of panel to be split
// into sub-panels, NS is the number of sub-panel on a panel and NF is the number of field cell to close to a particular
// panel.

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
