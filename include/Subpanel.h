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

//// Subpanel structure
// Structure containing information about sub-panels

#ifndef FPMV1_SUBPANEL_H
#define FPMV1_SUBPANEL_H

#include <vector>

struct Subpanel {

    int NS = 16; // number of sub-panel for each panel
    int NSs = 4; // chord/span-wise number of sub-panel for each panel
    double NC = 0.05; // sub-panel distance ratio in normal direction
    double LC = 0.75; // sub-panel distance ratio in longitudinal direction
    std::vector <int> sI; // surface panel indices
    std::vector <std::vector <int> > fI; // field panel indices
};

#endif //FPMV1_SUBPANEL_H
