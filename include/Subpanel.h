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
