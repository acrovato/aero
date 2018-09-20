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

#ifndef FPMV1_COMPUTE_SVARS_H
#define FPMV1_COMPUTE_SVARS_H

#include "Network.h"

void compute_sVars(bool symY, double sRef, double alpha, double Minf, Eigen::Vector3d &vInf,
                Eigen::MatrixX3d &vSigma, Network &bPan, double &cL, double &cD);

#endif //FPMV1_COMPUTE_SQ_H
