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

#ifndef FPMV1_INFCBB_H
#define FPMV1_INFCBB_H

using namespace std;

array<double,2> infcB(bool wakeFlag,int idTgt, int idSrc, double x, double y, double z,
                      double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4);

#endif //FPMV1_INFCBB_H
