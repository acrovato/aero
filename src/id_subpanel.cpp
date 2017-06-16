//// Subpanel identification
// Identify which field points are to close to which panel and prepare AIC matrices.
//
// I/O:
// - bPan: body panels (structure)
// - fPan: field panel (structure)
// - sp: subpanels (structure)
// - spAIC: subpanels AIC (structure)
//
// CONSIDER USING vector of structs instead of struct of vectors
// CONSIDER USING fixed size of boolean vector instead of spP
// CONSIDER USING octree (tree) to find points instead of double loop

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include "id_subpanel.h"

#define NDIM 3

using namespace std;
using namespace Eigen;

void id_subpanel(Network &bPan, Field &fPan, Subpanel &sp, Subpanel_AIC &spAIC) {

    // Temporary variables
    bool FLAG; // reset flag
    int ii; // index
    MatrixXd pV; // current panel vertices
    unsigned long size; // vector size
    pV.resize(4, NDIM);
    RowVector3d pC; // current field center
    Vector4d eL, d, dC, h, hC; // edge length, distance to edge & cut-off, distance to mid-edge & cut-off
    RowVector3d r1, r2, eXL, v; // distance to vertices, unit longitudinal edge vector, field point to mid-edge vector
    vector<int> vTemp; // temporary vector
    MatrixXd mTemp; // temporary matrix

    //// Begin
    cout << "Identifying sub-panels... " << flush;

    //// Identification of nearfield points
    for (int j = 0; j < bPan.nP; j++) {
        // Vertices coordinates & edges length & cut-off distance
        pV.row(0) = bPan.v0.row(j);
        pV.row(1) = bPan.v1.row(j);
        pV.row(2) = bPan.v2.row(j);
        pV.row(3) = bPan.v3.row(j);
        for (int k = 1; k <= 4; ++k) {
            int l = k % 4 + 1;
            eL(k-1) = (pV.row(l-1)-pV.row(k-1)).norm();
        }
        dC = sp.NC * eL;
        hC = sp.LC * eL;
        // Loop reset
        FLAG = 0;
        vTemp = {};
        for (int i = 0; i < fPan.nE; ++i) {
            ii = fPan.eIdx(i);
            // Target coordinates
            pC = fPan.CG.row(ii);
            // Distance to edges and mid-edges
            for (int k = 1; k <= 4; ++k) {
                int l = k % 4 + 1;
                r1 = pC - pV.row(k-1);
                r2 = pC - pV.row(l-1);
                eXL = (r2 - r1) / (r1 - r2).norm();
                v = (r1 + r2)/2;
                h(k-1) = abs(eXL.dot(v));
                d(k-1) = (r1.cross(r2)).norm() / (r1 - r2).norm();
            }
            // if target point is closer to edge than tolerance in perpendicular and longitudinal directions
            if (((d(0)-dC(0) <= 0) && (h(0)-hC(0) <= 0)) || ((d(1)-dC(1) <= 0) && (h(1)-hC(1) <= 0))
                || ((d(2)-dC(2) <= 0) && (h(2)-hC(2) <= 0)) || ((d(3)-dC(3) <= 0) && (h(3)-hC(3) <= 0))) {
                if (!FLAG) {
                    sp.sI.push_back(j);
                    FLAG = 1;
                }
                vTemp.push_back(ii);
            }
            if (i == fPan.nE-1 && FLAG)
                sp.fI.push_back(vTemp);
        }
    }

    //// Initialization of AIC matrices to 0
    for(int i = 0; i < sp.sI.size(); i++) {
        size = sp.fI[i].size();
        mTemp = MatrixXd::Zero(size, sp.NS);
        // Field
        spAIC.A.push_back(mTemp);
        spAIC.B.push_back(mTemp);
    }

    //// Control display
    cout << "Done!" << endl;
    cout << "Number of panels to split: " << sp.sI.size() << endl;
    cout << "Number of sub-panels per panel: " << sp.NS << endl;
    //for(int i = 0; i < sp.sI.size(); i++) {
    //    cout << sp.sI[i] << ": ";
    //    for(int j = 0; j < sp.fI[i].size(); j++) {
    //        cout << sp.fI[i][j] << ',';
    //    }
    //    cout << endl;
    //}
    cout << endl;
}