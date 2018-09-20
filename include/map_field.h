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

#ifndef FPMV1_MAP_FIELD_H
#define FPMV1_MAP_FIELD_H

#include "Field.h"
#include "Network.h"
#include "Numerical_CST.h"

void map_field(Eigen::MatrixX3d &sGrid, Numerical_CST &numC, Network &bPan, Field &fPan);

#endif //FPMV1_MAP_FIELD_H
