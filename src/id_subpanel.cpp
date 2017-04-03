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
#define NS 16
#define C 0.01

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
    Vector4d eL, d, dC; // edge length, distance to edge, cut-off
    RowVector3d r1, r2; // distance to vertices
    vector<int> vTemp; // temporary vector
    MatrixXd mTemp; // temporary matrix

    //// Begin
    cout << "Creating wake panels... ";

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
        dC = C * eL;
        // Loop reset
        FLAG = 0;
        vTemp = {};
        for (int i = 0; i < fPan.nE; ++i) {
            ii = fPan.eIdx(i);
            // Target coordinates
            pC = fPan.CG.row(ii);
            // Distance to edges
            for (int k = 1; k <= 4; ++k) {
                int l = k % 4 + 1;
                r1 = pC - pV.row(k-1);
                r2 = pC - pV.row(l-1);
                d(k-1) = (r1.cross(r2)).norm() / (r1 - r2).norm();
            }
            // if target point is closer to edge than ALPHA * (edge length)
            if ((d-dC).minCoeff() <= 0) {
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
        mTemp = MatrixXd::Zero(size, NS);
        // Field
        spAIC.Au.push_back(mTemp);
        spAIC.Av.push_back(mTemp);
        spAIC.Aw.push_back(mTemp);
        spAIC.Bu.push_back(mTemp);
        spAIC.Bv.push_back(mTemp);
        spAIC.Bw.push_back(mTemp);
        // Minigrid
        spAIC.AuXbwd.push_back(mTemp);
        spAIC.AvXbwd.push_back(mTemp);
        spAIC.AwXbwd.push_back(mTemp);
        spAIC.BuXbwd.push_back(mTemp);
        spAIC.BvXbwd.push_back(mTemp);
        spAIC.BwXbwd.push_back(mTemp);
        spAIC.AuXfwd.push_back(mTemp);
        spAIC.AvXfwd.push_back(mTemp);
        spAIC.AwXfwd.push_back(mTemp);
        spAIC.BuXfwd.push_back(mTemp);
        spAIC.BvXfwd.push_back(mTemp);
        spAIC.BwXfwd.push_back(mTemp);
        spAIC.AuYbwd.push_back(mTemp);
        spAIC.AvYbwd.push_back(mTemp);
        spAIC.AwYbwd.push_back(mTemp);
        spAIC.BuYbwd.push_back(mTemp);
        spAIC.BvYbwd.push_back(mTemp);
        spAIC.BwYbwd.push_back(mTemp);
        spAIC.AuYfwd.push_back(mTemp);
        spAIC.AvYfwd.push_back(mTemp);
        spAIC.AwYfwd.push_back(mTemp);
        spAIC.BuYfwd.push_back(mTemp);
        spAIC.BvYfwd.push_back(mTemp);
        spAIC.BwYfwd.push_back(mTemp);
        spAIC.AuZbwd.push_back(mTemp);
        spAIC.AvZbwd.push_back(mTemp);
        spAIC.AwZbwd.push_back(mTemp);
        spAIC.BuZbwd.push_back(mTemp);
        spAIC.BvZbwd.push_back(mTemp);
        spAIC.BwZbwd.push_back(mTemp);
        spAIC.AuZfwd.push_back(mTemp);
        spAIC.AvZfwd.push_back(mTemp);
        spAIC.AwZfwd.push_back(mTemp);
        spAIC.BuZfwd.push_back(mTemp);
        spAIC.BvZfwd.push_back(mTemp);
        spAIC.BwZfwd.push_back(mTemp);
    }

    //// Control display
    cout << "Done!" << endl;
    cout << "Number of panels to split: " << sp.sI.size() << endl;
    //for(int i = 0; i < sp.sI.size(); i++) {
    //    cout << sp.sI[i] << ", ";
    //    for(int j = 0; j < sp.fI[i].size(); j++) {
    //        cout << sp.fI[i][j] << ':';
    //    }
    //    cout << endl;
    //}
    cout << endl;
}