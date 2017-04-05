#ifndef FPMV1_SPLIT_PANEL_H
#define FPMV1_SPLIT_PANEL_H

std::array<Eigen::RowVectorXd,6> split_panel(int NS, int NSs, double x1, double x2, double x3, double x4,
                                 double y1, double y2, double y3, double y4,
                                 double x, double y, double z,
                                 double l0, double l1, double l2,
                                 double p0, double p1, double p2,
                                 double n0, double n1, double n2);

#endif //FPMV1_SPLIT_PANEL_H
