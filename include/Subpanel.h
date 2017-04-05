//// Subpanel
// Structure containing information about subpanels
//
// - NS: number of sub-panel on each panel
// - NS: number of sub-panel on each panel (along chord and span)
// - sI: vector of surface panel indices
// - fI: vector of vector of field indices (first dimension correspond to sI, second contains field indices)

#ifndef FPMV1_SUBPANEL_H
#define FPMV1_SUBPANEL_H

#include <vector>

struct Subpanel {

    const int NS = 16; // number of sub-panel for each panel //TODO: put this input file
    const int NSs = (const int) sqrt(NS); // chord/span-wise number of sub-panel for each panel
    std::vector <int> sI; // surface panel indices
    std::vector <std::vector <int> > fI; // field panel indices
};

#endif //FPMV1_SUBPANEL_H
