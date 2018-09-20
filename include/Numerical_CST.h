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

//// Numerical constants structure
// Structure containing numerical constants

#ifndef FPMV1_NUMERICAL_CST_H
#define FPMV1_NUMERICAL_CST_H

struct Numerical_CST {
    // Pre-processing
    double TOLB = 1e-3; // geometric tolerance on box (enclosing the geometry) size
    // Solver
    double RRED = 5; // order of magnitude of residual (relative change in field source) reduction
};

#endif //FPMV1_NUMERICAL_CST_H
