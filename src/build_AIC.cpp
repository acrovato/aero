//// AIC matrices building
// Build AIC matrices for body to body, body to field, field to field, and field to body interactions
//
// I/O:
// - symY: defines symmetry about Y axis
// - bPan: (network of) body panels (structure)
// - wPan: (network of) wake panels (structure)
// - fPan: field panels (structure)
// - b2bAIC: body to body AICs (structure)
// - b2fAIC: body to field AICs (structure)
// - f2fAIC: field to field AICs (structure)
// - f2bAIC: field to body AICs (structure)
// - sp: sub-panels (structure)
// - spAIC: body to field sub-panel AICs (structure)

#include <iostream>
#include <Eigen/Dense>
#include <array>
#include "build_AIC.h"
#include "infcB.h"
#include "infcFF.h"
#include "infcFB.h"
#include "split_panel.h"

using namespace std;
using namespace Eigen;

void build_AIC(bool symY, Network &bPan, Network &wPan, Field &fPan,
               Body_AIC &b2bAIC, Body_AIC &b2fAIC, Field2field_AIC &f2fAIC, Field2body_AIC &f2bAIC,
               Subpanel &sp, Subpanel_AIC &spAIC) {

    //// Initialization
    // Temporay variables
    int ii = 0, jj = 0; // counter
    bool FLAG = 0; // flag
    Vector3d dist; // Distance between panels
    Vector3d trsfColloc, trsfCorner1, trsfCorner2, trsfCorner3, trsfCorner4; // colloc (i) and corner (j) in local (j) coordinates

    Matrix3d glob2loc; // transformation matrix
    array<double, 2> coeffB; // body influence container

    array<RowVectorXd, 2> coeffSP; // body to field influence container for subpanel
    for (int i = 0; i < 2; i++)
        coeffSP[i].resize(sp.NS);

    array<double,2> xC, yC, zC; // cell vertices
    double x, y, z; // cell center
    double coeffFF; // field to field influence container
    array<double,3> coeffFB; // field to body influence container

    //// Begin
    cout << "Building AIC matrices... " << flush;

    // Resize AIC matrices (regular)
    b2bAIC.A.resize(bPan.nP, bPan.nP);
    b2bAIC.B.resize(bPan.nP, bPan.nP);
    b2fAIC.A.resize(fPan.nF, bPan.nP);
    b2fAIC.B.resize(fPan.nF, bPan.nP);
    f2fAIC.C.resize(fPan.nF, fPan.nF);
    f2bAIC.Cu.resize(bPan.nP, fPan.nF);
    f2bAIC.Cv.resize(bPan.nP, fPan.nF);
    f2bAIC.Cw.resize(bPan.nP, fPan.nF);

    //// Panel AIC
    for (int j = 0; j < bPan.nP; ++j) {
        glob2loc.row(0) = bPan.l.row(j);
        glob2loc.row(1) = bPan.p.row(j);
        glob2loc.row(2) = bPan.n.row(j);
        trsfCorner1 = bPan.v0.row(j).transpose();
        trsfCorner1 = glob2loc * trsfCorner1;
        trsfCorner2 = bPan.v1.row(j).transpose();
        trsfCorner2 = glob2loc * trsfCorner2;
        trsfCorner3 = bPan.v2.row(j).transpose();
        trsfCorner3 = glob2loc * trsfCorner3;
        trsfCorner4 = bPan.v3.row(j).transpose();
        trsfCorner4 = glob2loc * trsfCorner4;

        // panel to panel
        for (int i = 0; i < bPan.nP; ++i) {
            dist(0) = bPan.CG(i,0) - bPan.CG(j,0);
            dist(1) = bPan.CG(i,1) - bPan.CG(j,1);
            dist(2) = bPan.CG(i,2) - bPan.CG(j,2);
            dist = glob2loc * dist;
            trsfColloc = bPan.CG.row(i).transpose();
            trsfColloc = glob2loc * trsfColloc;

            coeffB = infcB(0,i,j, trsfColloc(0),trsfColloc(1),dist(2),trsfCorner1(0),trsfCorner1(1),
                           trsfCorner2(0),trsfCorner2(1),trsfCorner3(0),trsfCorner3(1),trsfCorner4(0),trsfCorner4(1));

            b2bAIC.A(i,j) = coeffB[0];
            b2bAIC.B(i,j) = coeffB[1];
        }

        // panel to field
        ii = 0; // loop reset
        FLAG = 0;
        for (int i = 0; i < fPan.nF; ++i) {
            // Split panel and store AIC
            if (sp.sI.size() != 0 && j == sp.sI[jj] && i == sp.fI[jj][ii]) {

                coeffSP = split_panel(j, i, trsfCorner1(0), trsfCorner2(0), trsfCorner3(0), trsfCorner4(0),
                                      trsfCorner1(1), trsfCorner2(1), trsfCorner3(1), trsfCorner4(1),
                                      bPan, fPan, sp);
                spAIC.A[jj].row(ii) = coeffSP[0];
                spAIC.B[jj].row(ii) = coeffSP[1];
                // Set global AIC to 0
                b2fAIC.A(i,j) = 0;
                b2fAIC.B(i,j) = 0;

                // Set flags
                if (ii == 0)
                    FLAG = 1;
                if (ii < sp.fI[jj].size()-1)
                    ii++;
            }
            else {
                dist(0) = fPan.CG(i,0) - bPan.CG(j,0);
                dist(1) = fPan.CG(i,1) - bPan.CG(j,1);
                dist(2) = fPan.CG(i,2) - bPan.CG(j,2);
                dist = glob2loc * dist;
                trsfColloc = fPan.CG.row(i).transpose();
                trsfColloc = glob2loc * trsfColloc;

                coeffB = infcB(0, 0, 1, trsfColloc(0), trsfColloc(1), dist(2), trsfCorner1(0), trsfCorner1(1),
                                 trsfCorner2(0), trsfCorner2(1), trsfCorner3(0), trsfCorner3(1), trsfCorner4(0), trsfCorner4(1));

                b2fAIC.A(i,j) = coeffB[0];
                b2fAIC.B(i,j) = coeffB[1];
            }
            if (i == fPan.nF-1 && FLAG && jj < sp.sI.size()-1)
                jj++;
        }
    }

    // Panel symmetry
    if (symY) {
        for (int j = 0; j < bPan.nP; ++j) {
            glob2loc.row(0) = bPan.l.row(j);
            glob2loc.row(1) = -bPan.p.row(j);
            glob2loc.row(2) = bPan.n.row(j);
            glob2loc.col(1) = -glob2loc.col(1);
            trsfCorner1 = bPan.v0.row(j).transpose();
            trsfCorner1(1) = -trsfCorner1(1);
            trsfCorner1 = glob2loc * trsfCorner1;
            trsfCorner2 = bPan.v1.row(j).transpose();
            trsfCorner2(1) = -trsfCorner2(1);
            trsfCorner2 = glob2loc * trsfCorner2;
            trsfCorner3 = bPan.v2.row(j).transpose();
            trsfCorner3(1) = -trsfCorner3(1);
            trsfCorner3 = glob2loc * trsfCorner3;
            trsfCorner4 = bPan.v3.row(j).transpose();
            trsfCorner4(1) = -trsfCorner4(1);
            trsfCorner4 = glob2loc * trsfCorner4;

            // panel to panel
            for (int i = 0; i < bPan.nP; ++i) {
                dist(0) = bPan.CG(i,0) - bPan.CG(j,0);
                dist(1) = bPan.CG(i,1) + bPan.CG(j,1);
                dist(2) = bPan.CG(i,2) - bPan.CG(j,2);
                dist = glob2loc * dist;
                trsfColloc = bPan.CG.row(i).transpose();
                trsfColloc = glob2loc * trsfColloc;

                coeffB = infcB(0,0,1,trsfColloc(0),trsfColloc(1),dist(2),trsfCorner1(0),trsfCorner1(1),
                               trsfCorner2(0),trsfCorner2(1),trsfCorner3(0),trsfCorner3(1),trsfCorner4(0),trsfCorner4(1));

                b2bAIC.A(i,j) -= coeffB[0]; // minus sign because corner order is trigonometric!
                b2bAIC.B(i,j) -= coeffB[1];
            }

            // panel to field
            for (int i = 0; i < fPan.nF; ++i) {
                dist(0) = fPan.CG(i, 0) - bPan.CG(j, 0);
                dist(1) = fPan.CG(i, 1) + bPan.CG(j, 1);
                dist(2) = fPan.CG(i, 2) - bPan.CG(j, 2);
                dist = glob2loc * dist;
                trsfColloc = fPan.CG.row(i).transpose();
                trsfColloc = glob2loc * trsfColloc;

                coeffB = infcB(0, 0, 1, trsfColloc(0), trsfColloc(1), dist(2), trsfCorner1(0), trsfCorner1(1),
                                trsfCorner2(0), trsfCorner2(1), trsfCorner3(0), trsfCorner3(1), trsfCorner4(0), trsfCorner4(1));

                b2fAIC.A(i,j) -= coeffB[0];
                b2fAIC.B(i,j) -= coeffB[1];
            }
        }
    }

    //// Wake panel AIC
    // each panel i_th is influenced by each wake panel (total number m) whose mu_M = mu_1M - mu_NM
    // on each row i, columns 0 and n-1 of each wake are modified (j*n and (j+1)*n -1)
    for (int j = 0; j < wPan.nS_; ++j) {
        glob2loc.row(0) = wPan.l.row(j);
        glob2loc.row(1) = wPan.p.row(j);
        glob2loc.row(2) = wPan.n.row(j);
        trsfCorner1 = wPan.v0.row(j).transpose();
        trsfCorner1 = glob2loc * trsfCorner1;
        trsfCorner2 = wPan.v1.row(j).transpose();
        trsfCorner2 = glob2loc * trsfCorner2;
        trsfCorner3 = wPan.v2.row(j).transpose();
        trsfCorner3 = glob2loc * trsfCorner3;
        trsfCorner4 = wPan.v3.row(j).transpose();
        trsfCorner4 = glob2loc * trsfCorner4;

        // panel to panel
        for (int i = 0; i < bPan.nP; ++i) {
            dist(0) = bPan.CG(i, 0) - wPan.CG(j, 0);
            dist(1) = bPan.CG(i, 1) - wPan.CG(j, 1);
            dist(2) = bPan.CG(i, 2) - wPan.CG(j, 2);
            dist = glob2loc * dist;
            trsfColloc = bPan.CG.row(i).transpose();
            trsfColloc = glob2loc * trsfColloc;

            coeffB = infcB(1,i,j,trsfColloc(0),trsfColloc(1),dist(2),trsfCorner1(0),trsfCorner1(1),
                           trsfCorner2(0),trsfCorner2(1),trsfCorner3(0),trsfCorner3(1),trsfCorner4(0),trsfCorner4(1));

            b2bAIC.A(i, j * bPan.nC_) += coeffB[0];
            b2bAIC.A(i, (j + 1) * bPan.nC_ - 1) -= coeffB[0];
        }

        // panel to field
        for (int i = 0; i < fPan.nF; ++i) {
            dist(0) = fPan.CG(i, 0) - wPan.CG(j, 0);
            dist(1) = fPan.CG(i, 1) - wPan.CG(j, 1);
            dist(2) = fPan.CG(i, 2) - wPan.CG(j, 2);
            dist = glob2loc * dist;
            trsfColloc = fPan.CG.row(i).transpose();
            trsfColloc = glob2loc * trsfColloc;

            coeffB = infcB(1, 0, 1, trsfColloc(0), trsfColloc(1), dist(2), trsfCorner1(0), trsfCorner1(1),
                            trsfCorner2(0), trsfCorner2(1), trsfCorner3(0), trsfCorner3(1), trsfCorner4(0), trsfCorner4(1));

            b2fAIC.A(i, j * bPan.nC_) += coeffB[0];
            b2fAIC.A(i, (j + 1) * bPan.nC_ - 1) -= coeffB[0];
        }
    }

    // Wake symmetry
    if (symY) {
        for (int j = 0; j < wPan.nS_; ++j) {
            glob2loc.row(0) = wPan.l.row(j);
            glob2loc.row(1) = -wPan.p.row(j);
            glob2loc.row(2) = wPan.n.row(j);
            glob2loc.col(1) = -glob2loc.col(1);
            trsfCorner1 = wPan.v0.row(j).transpose();
            trsfCorner1(1) = -trsfCorner1(1);
            trsfCorner1 = glob2loc * trsfCorner1;
            trsfCorner2 = wPan.v1.row(j).transpose();
            trsfCorner2(1) = -trsfCorner2(1);
            trsfCorner2 = glob2loc * trsfCorner2;
            trsfCorner3 = wPan.v2.row(j).transpose();
            trsfCorner3(1) = -trsfCorner3(1);
            trsfCorner3 = glob2loc * trsfCorner3;
            trsfCorner4 = wPan.v3.row(j).transpose();
            trsfCorner4(1) = -trsfCorner4(1);
            trsfCorner4 = glob2loc * trsfCorner4;

            // panel to panel
            for (int i = 0; i < bPan.nP; ++i) {
                dist(0) = bPan.CG(i, 0) - wPan.CG(j, 0);
                dist(1) = bPan.CG(i, 1) + wPan.CG(j, 1);
                dist(2) = bPan.CG(i, 2) - wPan.CG(j, 2);
                dist = glob2loc * dist;
                trsfColloc = bPan.CG.row(i).transpose();
                trsfColloc = glob2loc * trsfColloc;

                coeffB = infcB(1,0,1,trsfColloc(0),trsfColloc(1),dist(2),trsfCorner1(0),trsfCorner1(1),
                               trsfCorner2(0),trsfCorner2(1),trsfCorner3(0),trsfCorner3(1),trsfCorner4(0),trsfCorner4(1));

                b2bAIC.A(i, j * bPan.nC_) -= coeffB[0];
                b2bAIC.A(i, (j + 1) * bPan.nC_ - 1) += coeffB[0];
            }

            // panel to field
            for (int i = 0; i < fPan.nF; ++i) {
                dist(0) = fPan.CG(i, 0) - wPan.CG(j, 0);
                dist(1) = fPan.CG(i, 1) + wPan.CG(j, 1);
                dist(2) = fPan.CG(i, 2) - wPan.CG(j, 2);
                dist = glob2loc * dist;
                trsfColloc = fPan.CG.row(i).transpose();
                trsfColloc = glob2loc * trsfColloc;

                coeffB = infcB(1, 0, 1, trsfColloc(0), trsfColloc(1), dist(2), trsfCorner1(0), trsfCorner1(1),
                                trsfCorner2(0), trsfCorner2(1), trsfCorner3(0), trsfCorner3(1), trsfCorner4(0), trsfCorner4(1));

                b2fAIC.A(i, j * bPan.nC_) -= coeffB[0];
                b2fAIC.A(i, (j + 1) * bPan.nC_ - 1) += coeffB[0];
            }
        }
    }

    //// Field AIC
    for (int j = 0; j < fPan.nF; ++j) {
        xC[0] = fPan.vX(j,0) - fPan.CG(j,0);
        xC[1] = fPan.vX(j,1) - fPan.CG(j,0);
        yC[0] = fPan.vY(j,0) - fPan.CG(j,1);
        yC[1] = fPan.vY(j,1) - fPan.CG(j,1);
        zC[0] = fPan.vZ(j,0) - fPan.CG(j,2);
        zC[1] = fPan.vZ(j,1) - fPan.CG(j,2);

        // Field to field (diagonal)
        x = fPan.CG(j,0) - fPan.CG(j,0);
        y = fPan.CG(j,1) - fPan.CG(j,1);
        z = fPan.CG(j,2) - fPan.CG(j,2);
        coeffFF = infcFF(x, y, z, xC, yC, zC);
        f2fAIC.C(j,j) = coeffFF;
        // Field to field
        for (int i = j+1; i < fPan.nF; ++i) {
            x = fPan.CG(i,0) - fPan.CG(j,0);
            y = fPan.CG(i,1) - fPan.CG(j,1);
            z = fPan.CG(i,2) - fPan.CG(j,2);
            coeffFF = infcFF(x, y, z, xC, yC, zC);
            f2fAIC.C(i,j) = coeffFF;
            f2fAIC.C(j,i) = coeffFF;
        }

        // Field to panel
        for (int i = 0; i < bPan.nP; ++i) {
            x = bPan.CG(i,0) - fPan.CG(j,0);
            y = bPan.CG(i,1) - fPan.CG(j,1);
            z = bPan.CG(i,2) - fPan.CG(j,2);
            coeffFB = infcFB(x, y, z, xC, yC, zC);
            f2bAIC.Cu(i,j) = coeffFB[0];
            f2bAIC.Cv(i,j) = coeffFB[1];
            f2bAIC.Cw(i,j) = coeffFB[2];
        }
    }

    // Field symmetry
    if (symY) {
        for (int j = 0; j < fPan.nF; ++j) {
            xC[0] = fPan.vX(j,0) - fPan.CG(j,0);
            xC[1] = fPan.vX(j,1) - fPan.CG(j,0);
            yC[1] = -fPan.vY(j,0) + fPan.CG(j,1); // Y corner order inverted to preserve sign
            yC[0] = -fPan.vY(j,1) + fPan.CG(j,1);
            zC[0] = fPan.vZ(j,0) - fPan.CG(j,2);
            zC[1] = fPan.vZ(j,1) - fPan.CG(j,2);

            // Field to field
            for (int i = 0; i < fPan.nF; ++i) {
                x = fPan.CG(i,0) - fPan.CG(j,0);
                y = fPan.CG(i,1) + fPan.CG(j,1);
                z = fPan.CG(i,2) - fPan.CG(j,2);
                coeffFF = infcFF(x, y, z, xC, yC, zC);
                f2fAIC.C(i,j) += coeffFF;
            }

            // Field to panel
            for (int i = 0; i < bPan.nP; ++i) {
                x = bPan.CG(i,0) - fPan.CG(j,0);
                y = bPan.CG(i,1) + fPan.CG(j,1);
                z = bPan.CG(i,2) - fPan.CG(j,2);
                coeffFB = infcFB(x, y, z, xC, yC, zC);
                f2bAIC.Cu(i,j) += coeffFB[0];
                f2bAIC.Cv(i,j) += coeffFB[1];
                f2bAIC.Cw(i,j) += coeffFB[2];
            }
        }
    }

    //// Control display
    cout << "Done!" << endl << endl;
}