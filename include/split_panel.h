#ifndef FPMV1_SPLIT_PANEL_H
#define FPMV1_SPLIT_PANEL_H

#include "Network.h"
#include "Field.h"
#include "Subpanel.h"

std::array<Eigen::RowVectorXd,6> split_panel(int idP, int idF,
                                              double x0, double x1, double x2, double x3, double y0, double y1, double y2, double y3,
                                              Network &bPan, Field &fPan, Subpanel &sp);

#endif //FPMV1_SPLIT_PANEL_H
