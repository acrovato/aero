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

//// Field to body AIC structure
// Structure containing field to body AIC matrices

#ifndef FPMV1_FIELD2BODY_AIC_H
#define FPMV1_FIELD2BODY_AIC_H

#include <Eigen/Dense>

struct Field2body_AIC {
    Eigen::MatrixXd Cu; // source matrix for x-component of velocity
    Eigen::MatrixXd Cv; // source matrix for y-component of velocity
    Eigen::MatrixXd Cw; // source matrix for z-component of velocity
};

#endif //FPMV1_FIELD2BODY_AIC_H
