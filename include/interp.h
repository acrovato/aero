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

#ifndef FPMV1_INTERP_H
#define FPMV1_INTERP_H

double interp(double x0, double y0, double z0, double x1, double y1, double z1,
              double x2, double y2, double z2, double x3, double y3, double z3,
              double s0, double s1, double s2, double s3,
              double x, double y, double z);

#endif //FPMV1_INTERP_H
