
/* Copyright (C) 2018 Adrien Crovato */

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
