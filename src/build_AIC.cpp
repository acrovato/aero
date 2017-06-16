//// Build AIC
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
// - mgAIC: body to field AICs for minigrid (structure)

#include <iostream>
#include <Eigen/Dense>
#include <array>
#include "build_AIC.h"
#include "infcBB.h"
#include "infcBF.h"
#include "infcF.h"
#include "split_panel.h"

using namespace std;
using namespace Eigen;

void build_AIC(bool symY, Network &bPan, Network &wPan, Field &fPan,
               Body_AIC &b2bAIC, Body2field_AIC &b2fAIC, Field_AIC &f2fAIC, Field_AIC &f2bAIC,
               Subpanel &sp, Subpanel_AIC &spAIC) {

    //// Initialization
    // Temporay variables
    int ii = 0, jj = 0; // counter
    bool FLAG = 0; // flag
    Vector3d dist; // Distance between panels
    Vector3d trsfColloc, trsfCorner1, trsfCorner2, trsfCorner3, trsfCorner4; // colloc (i) and corner (j) in local (j) coordinates

    Matrix3d glob2loc; // transformation matrix
    array<double, 2> coeffBB; // body to body influence container
    array<double, 6> coeffBF; // body to field influence container

    array<RowVectorXd, 6> coeffSP; // body to field influence container for subpanel
    for (int i = 0; i < 6; i++)
        coeffSP[i].resize(sp.NS);

    array<double,2> xC, yC, zC; // cell vertices
    double x, y, z; // cell center
    array<double,3> coeffF; // field influence container

    //// Begin
    cout << "Building AIC matrices... " << flush;

    // Resize AIC matrices (regular)
    b2bAIC.A.resize(bPan.nP, bPan.nP);
    b2bAIC.B.resize(bPan.nP, bPan.nP);
    b2fAIC.Au.resize(fPan.nF, bPan.nP);
    b2fAIC.Av.resize(fPan.nF, bPan.nP);
    b2fAIC.Aw.resize(fPan.nF, bPan.nP);
    b2fAIC.Bu.resize(fPan.nF, bPan.nP);
    b2fAIC.Bv.resize(fPan.nF, bPan.nP);
    b2fAIC.Bw.resize(fPan.nF, bPan.nP);
    f2fAIC.Cu.resize(fPan.nF, fPan.nF);
    f2fAIC.Cv.resize(fPan.nF, fPan.nF);
    f2fAIC.Cw.resize(fPan.nF, fPan.nF);
    f2bAIC.Cu.resize(bPan.nP, fPan.nF);
    f2bAIC.Cv.resize(bPan.nP, fPan.nF);
    f2bAIC.Cw.resize(bPan.nP, fPan.nF);

    // TODO: rework build_AIC. Group parts with minigrid.

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

            coeffBB = infcBB(0,i,j, trsfColloc(0),trsfColloc(1),dist(2),trsfCorner1(0),trsfCorner1(1),
                           trsfCorner2(0),trsfCorner2(1),trsfCorner3(0),trsfCorner3(1),trsfCorner4(0),trsfCorner4(1));

            b2bAIC.A(i,j) = coeffBB[0];
            b2bAIC.B(i,j) = coeffBB[1];
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
                /// Center
                spAIC.Au[jj].row(ii) = coeffSP[0];
                spAIC.Av[jj].row(ii) = coeffSP[1];
                spAIC.Aw[jj].row(ii) = coeffSP[2];
                spAIC.Bu[jj].row(ii) = coeffSP[3];
                spAIC.Bv[jj].row(ii) = coeffSP[4];
                spAIC.Bw[jj].row(ii) = coeffSP[5];
                // Set global AIC to 0
                b2fAIC.Au(i,j) = 0;
                b2fAIC.Av(i,j) = 0;
                b2fAIC.Aw(i,j) = 0;
                b2fAIC.Bu(i,j) = 0;
                b2fAIC.Bv(i,j) = 0;
                b2fAIC.Bw(i,j) = 0;

                // Set flags
                if (ii == 0)
                    FLAG = 1;
                if (ii < sp.fI[jj].size()-1)
                    ii++;
            }
            else {
                // Center
                dist(0) = fPan.CG(i,0) - bPan.CG(j,0);
                dist(1) = fPan.CG(i,1) - bPan.CG(j,1);
                dist(2) = fPan.CG(i,2) - bPan.CG(j,2);
                dist = glob2loc * dist;
                trsfColloc = fPan.CG.row(i).transpose();
                trsfColloc = glob2loc * trsfColloc;

                coeffBF = infcBF(0, trsfCorner1(0), trsfCorner2(0), trsfCorner3(0), trsfCorner4(0),
                                 trsfCorner1(1), trsfCorner2(1), trsfCorner3(1), trsfCorner4(1),
                                 trsfColloc(0), trsfColloc(1), dist(2),
                                 bPan.l(j, 0), bPan.l(j, 1), bPan.l(j, 2), bPan.p(j, 0), bPan.p(j, 1), bPan.p(j, 2),
                                 bPan.n(j, 0), bPan.n(j, 1), bPan.n(j, 2));
                b2fAIC.Au(i,j) = coeffBF[0];
                b2fAIC.Av(i,j) = coeffBF[1];
                b2fAIC.Aw(i,j) = coeffBF[2];
                b2fAIC.Bu(i,j) = coeffBF[3];
                b2fAIC.Bv(i,j) = coeffBF[4];
                b2fAIC.Bw(i,j) = coeffBF[5];
            }
            if (i == fPan.nF-1 && FLAG && jj < sp.sI.size()-1)
                jj++;
        }
    }

    // TODO check panel vectors symmetry for complex wings
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

                coeffBB = infcBB(0,0,1,trsfColloc(0),trsfColloc(1),dist(2),trsfCorner1(0),trsfCorner1(1),
                               trsfCorner2(0),trsfCorner2(1),trsfCorner3(0),trsfCorner3(1),trsfCorner4(0),trsfCorner4(1));

                b2bAIC.A(i,j) -= coeffBB[0]; // minus sign because corner order is trigonometric!
                b2bAIC.B(i,j) -= coeffBB[1];
            }

            // panel to field
            for (int i = 0; i < fPan.nF; ++i) {
                dist(0) = fPan.CG(i, 0) - bPan.CG(j, 0);
                dist(1) = fPan.CG(i, 1) + bPan.CG(j, 1);
                dist(2) = fPan.CG(i, 2) - bPan.CG(j, 2);
                dist = glob2loc * dist;
                trsfColloc = fPan.CG.row(i).transpose();
                trsfColloc = glob2loc * trsfColloc;

                coeffBF = infcBF(0, trsfCorner1(0) ,trsfCorner2(0), trsfCorner3(0), trsfCorner4(0),
                                 trsfCorner1(1), trsfCorner2(1), trsfCorner3(1), trsfCorner4(1),
                                 trsfColloc(0), trsfColloc(1), dist(2),
                                 bPan.l(j,0), -bPan.l(j,1), bPan.l(j,2), -bPan.p(j,0), bPan.p(j,1), -bPan.p(j,2),
                                 bPan.n(j,0), -bPan.n(j,1), bPan.n(j,2));

                b2fAIC.Au(i,j) -= coeffBF[0];
                b2fAIC.Av(i,j) -= coeffBF[1];
                b2fAIC.Aw(i,j) -= coeffBF[2];
                b2fAIC.Bu(i,j) -= coeffBF[3];
                b2fAIC.Bv(i,j) -= coeffBF[4];
                b2fAIC.Bw(i,j) -= coeffBF[5];
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

            coeffBB = infcBB(1,i,j,trsfColloc(0),trsfColloc(1),dist(2),trsfCorner1(0),trsfCorner1(1),
                           trsfCorner2(0),trsfCorner2(1),trsfCorner3(0),trsfCorner3(1),trsfCorner4(0),trsfCorner4(1));

            b2bAIC.A(i, j * bPan.nC_) += coeffBB[0];
            b2bAIC.A(i, (j + 1) * bPan.nC_ - 1) -= coeffBB[0];
        }

        // panel to field
        for (int i = 0; i < fPan.nF; ++i) {
            dist(0) = fPan.CG(i, 0) - wPan.CG(j, 0);
            dist(1) = fPan.CG(i, 1) - wPan.CG(j, 1);
            dist(2) = fPan.CG(i, 2) - wPan.CG(j, 2);
            dist = glob2loc * dist;
            trsfColloc = fPan.CG.row(i).transpose();
            trsfColloc = glob2loc * trsfColloc;

            coeffBF = infcBF(1, trsfCorner1(0) ,trsfCorner2(0), trsfCorner3(0), trsfCorner4(0),
                             trsfCorner1(1), trsfCorner2(1), trsfCorner3(1), trsfCorner4(1),
                             trsfColloc(0), trsfColloc(1), dist(2),
                             wPan.l(j,0), wPan.l(j,1), wPan.l(j,2), wPan.p(j,0), wPan.p(j,1), wPan.p(j,2),
                             wPan.n(j,0), wPan.n(j,1), wPan.n(j,2));

            b2fAIC.Au(i, j * bPan.nC_) += coeffBF[0];
            b2fAIC.Au(i, (j + 1) * bPan.nC_ - 1) -= coeffBF[0];
            b2fAIC.Av(i, j * bPan.nC_) += coeffBF[1];
            b2fAIC.Av(i, (j + 1) * bPan.nC_ - 1) -= coeffBF[1];
            b2fAIC.Aw(i, j * bPan.nC_) += coeffBF[2];
            b2fAIC.Aw(i, (j + 1) * bPan.nC_ - 1) -= coeffBF[2];
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

                coeffBB = infcBB(1,0,1,trsfColloc(0),trsfColloc(1),dist(2),trsfCorner1(0),trsfCorner1(1),
                               trsfCorner2(0),trsfCorner2(1),trsfCorner3(0),trsfCorner3(1),trsfCorner4(0),trsfCorner4(1));

                b2bAIC.A(i, j * bPan.nC_) -= coeffBB[0];
                b2bAIC.A(i, (j + 1) * bPan.nC_ - 1) += coeffBB[0];
            }

            // panel to field
            for (int i = 0; i < fPan.nF; ++i) {
                dist(0) = fPan.CG(i, 0) - wPan.CG(j, 0);
                dist(1) = fPan.CG(i, 1) + wPan.CG(j, 1);
                dist(2) = fPan.CG(i, 2) - wPan.CG(j, 2);
                dist = glob2loc * dist;
                trsfColloc = fPan.CG.row(i).transpose();
                trsfColloc = glob2loc * trsfColloc;

                coeffBF = infcBF(1, trsfCorner1(0) ,trsfCorner2(0), trsfCorner3(0), trsfCorner4(0),
                                 trsfCorner1(1), trsfCorner2(1), trsfCorner3(1), trsfCorner4(1),
                                 trsfColloc(0), trsfColloc(1), dist(2),
                                 wPan.l(j,0), -wPan.l(j,1), wPan.l(j,2), -wPan.p(j,0), wPan.p(j,1), -wPan.p(j,2),
                                 wPan.n(j,0), -wPan.n(j,1), wPan.n(j,2));

                b2fAIC.Au(i, j * bPan.nC_) -= coeffBF[0];
                b2fAIC.Au(i, (j + 1) * bPan.nC_ - 1) += coeffBF[0];
                b2fAIC.Av(i, j * bPan.nC_) -= coeffBF[1];
                b2fAIC.Av(i, (j + 1) * bPan.nC_ - 1) += coeffBF[1];
                b2fAIC.Aw(i, j * bPan.nC_) -= coeffBF[2];
                b2fAIC.Aw(i, (j + 1) * bPan.nC_ - 1) += coeffBF[2];
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
        coeffF = infcF(x, y, z, xC, yC, zC);
        f2fAIC.Cu(j,j) = coeffF[0];
        f2fAIC.Cv(j,j) = coeffF[1];
        f2fAIC.Cw(j,j) = coeffF[2];
        // Field to field
        for (int i = j+1; i < fPan.nF; ++i) {
            x = fPan.CG(i,0) - fPan.CG(j,0);
            y = fPan.CG(i,1) - fPan.CG(j,1);
            z = fPan.CG(i,2) - fPan.CG(j,2);
            coeffF = infcF(x, y, z, xC, yC, zC);
            f2fAIC.Cu(i,j) = coeffF[0];
            f2fAIC.Cv(i,j) = coeffF[1];
            f2fAIC.Cw(i,j) = coeffF[2];
            f2fAIC.Cu(j,i) = -coeffF[0];
            f2fAIC.Cv(j,i) = -coeffF[1];
            f2fAIC.Cw(j,i) = -coeffF[2];
        }

        // Field to panel
        for (int i = 0; i < bPan.nP; ++i) {
            x = bPan.CG(i,0) - fPan.CG(j,0);
            y = bPan.CG(i,1) - fPan.CG(j,1);
            z = bPan.CG(i,2) - fPan.CG(j,2);
            coeffF = infcF(x, y, z, xC, yC, zC);
            f2bAIC.Cu(i,j) = coeffF[0];
            f2bAIC.Cv(i,j) = coeffF[1];
            f2bAIC.Cw(i,j) = coeffF[2];
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
                coeffF = infcF(x, y, z, xC, yC, zC);
                f2fAIC.Cu(i,j) += coeffF[0];
                f2fAIC.Cv(i,j) += coeffF[1];
                f2fAIC.Cw(i,j) += coeffF[2];
            }

            // Field to panel
            for (int i = 0; i < bPan.nP; ++i) {
                x = bPan.CG(i,0) - fPan.CG(j,0);
                y = bPan.CG(i,1) + fPan.CG(j,1);
                z = bPan.CG(i,2) - fPan.CG(j,2);
                coeffF = infcF(x, y, z, xC, yC, zC);
                f2bAIC.Cu(i,j) += coeffF[0];
                f2bAIC.Cv(i,j) += coeffF[1];
                f2bAIC.Cw(i,j) += coeffF[2];
            }
        }
    }

    //// Control display
    cout << "Done!" << endl << endl;
}