#ifndef FPMV1_AIC_H
#define FPMV1_AIC_H

#include "Network.h"
#include "Field.h"
#include "Body_AIC.h"
#include "Body2field_AIC.h"
#include "Field_AIC.h"
#include "Subpanel.h"
#include "Subpanel_AIC.h"

void build_AIC(bool symY, Network &bPan, Network &wPan, Field &fPan,
               Body_AIC &b2bAIC, Body2field_AIC &b2fAIC, Field_AIC &f2fAIC, Field_AIC &f2bAIC,
               Subpanel &sp, Subpanel_AIC &spAIC);

#endif //FPMV1_AIC_H
