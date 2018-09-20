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

#ifndef FPMV1_AIC_H
#define FPMV1_AIC_H

#include "Network.h"
#include "Field.h"
#include "Body_AIC.h"
#include "Field2field_AIC.h"
#include "Field2body_AIC.h"
#include "Subpanel.h"
#include "Subpanel_AIC.h"

void build_AIC(bool symY, Network &bPan, Network &wPan, Field &fPan,
               Body_AIC &b2bAIC, Body_AIC &b2fAIC, Field2field_AIC &f2fAIC, Field2body_AIC &f2bAIC,
               Subpanel &sp, Subpanel_AIC &spAIC);

#endif //FPMV1_AIC_H
