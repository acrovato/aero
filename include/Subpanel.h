//// Subpanel
// Structure containing information about subpanels
//
// - sI: vector of surface panel indices
// - fI: vector of vector of field indices (first dimension correspond to sI, second contains field indices)

#ifndef FPMV1_SUBPANEL_H
#define FPMV1_SUBPANEL_H

#include <vector>

struct Subpanel {

    std::vector <int> sI; // surface panel indices
    std::vector <std::vector <int> > fI; // field panel indices
};

#endif //FPMV1_SUBPANEL_H
