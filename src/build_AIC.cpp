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

#define NS 16

using namespace std;
using namespace Eigen;

void build_AIC(bool symY, Network &bPan, Network &wPan, Field &fPan,
               Body_AIC &b2bAIC, Body2field_AIC &b2fAIC, Field_AIC &f2fAIC, Field_AIC &f2bAIC, Minigrid_AIC &mgAIC,
               Subpanel &sp, Subpanel_AIC &spAIC) {

    // Temporay variables
    int ii = 0, jj = 0; // counter
    bool FLAG = 0; // flag
    Vector3d dist; // Distance between panels
    Vector3d trsfColloc, trsfCorner1, trsfCorner2, trsfCorner3, trsfCorner4; // colloc (i) and corner (j) in local (j) coordinates

    Matrix3d glob2loc; // transformation matrix
    array<double, 2> coeffBB; // body to body influence container
    array<double, 6> coeffBF; // body to field influence container

    array<RowVectorXd, 6> coeffSP; // body to field influence container for subpanel
    coeffSP[0].resize(NS);
    coeffSP[1].resize(NS);
    coeffSP[2].resize(NS);
    coeffSP[3].resize(NS);
    coeffSP[4].resize(NS);
    coeffSP[5].resize(NS);

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
    // Resize AIC matrices (minigrid)
    mgAIC.AuXbwd.resize(fPan.nF, bPan.nP);
    mgAIC.AvXbwd.resize(fPan.nF, bPan.nP);
    mgAIC.AwXbwd.resize(fPan.nF, bPan.nP);
    mgAIC.BuXbwd.resize(fPan.nF, bPan.nP);
    mgAIC.BvXbwd.resize(fPan.nF, bPan.nP);
    mgAIC.BwXbwd.resize(fPan.nF, bPan.nP);
    mgAIC.CuXbwd.resize(fPan.nF, fPan.nF);
    mgAIC.CvXbwd.resize(fPan.nF, fPan.nF);
    mgAIC.CwXbwd.resize(fPan.nF, fPan.nF);
    mgAIC.AuXfwd.resize(fPan.nF, bPan.nP);
    mgAIC.AvXfwd.resize(fPan.nF, bPan.nP);
    mgAIC.AwXfwd.resize(fPan.nF, bPan.nP);
    mgAIC.BuXfwd.resize(fPan.nF, bPan.nP);
    mgAIC.BvXfwd.resize(fPan.nF, bPan.nP);
    mgAIC.BwXfwd.resize(fPan.nF, bPan.nP);
    mgAIC.CuXfwd.resize(fPan.nF, fPan.nF);
    mgAIC.CvXfwd.resize(fPan.nF, fPan.nF);
    mgAIC.CwXfwd.resize(fPan.nF, fPan.nF);
    mgAIC.AuYbwd.resize(fPan.nF, bPan.nP);
    mgAIC.AvYbwd.resize(fPan.nF, bPan.nP);
    mgAIC.AwYbwd.resize(fPan.nF, bPan.nP);
    mgAIC.BuYbwd.resize(fPan.nF, bPan.nP);
    mgAIC.BvYbwd.resize(fPan.nF, bPan.nP);
    mgAIC.BwYbwd.resize(fPan.nF, bPan.nP);
    mgAIC.CuYbwd.resize(fPan.nF, fPan.nF);
    mgAIC.CvYbwd.resize(fPan.nF, fPan.nF);
    mgAIC.CwYbwd.resize(fPan.nF, fPan.nF);
    mgAIC.AuYfwd.resize(fPan.nF, bPan.nP);
    mgAIC.AvYfwd.resize(fPan.nF, bPan.nP);
    mgAIC.AwYfwd.resize(fPan.nF, bPan.nP);
    mgAIC.BuYfwd.resize(fPan.nF, bPan.nP);
    mgAIC.BvYfwd.resize(fPan.nF, bPan.nP);
    mgAIC.BwYfwd.resize(fPan.nF, bPan.nP);
    mgAIC.CuYfwd.resize(fPan.nF, fPan.nF);
    mgAIC.CvYfwd.resize(fPan.nF, fPan.nF);
    mgAIC.CwYfwd.resize(fPan.nF, fPan.nF);
    mgAIC.AuZbwd.resize(fPan.nF, bPan.nP);
    mgAIC.AvZbwd.resize(fPan.nF, bPan.nP);
    mgAIC.AwZbwd.resize(fPan.nF, bPan.nP);
    mgAIC.BuZbwd.resize(fPan.nF, bPan.nP);
    mgAIC.BvZbwd.resize(fPan.nF, bPan.nP);
    mgAIC.BwZbwd.resize(fPan.nF, bPan.nP);
    mgAIC.CuZbwd.resize(fPan.nF, fPan.nF);
    mgAIC.CvZbwd.resize(fPan.nF, fPan.nF);
    mgAIC.CwZbwd.resize(fPan.nF, fPan.nF);
    mgAIC.AuZfwd.resize(fPan.nF, bPan.nP);
    mgAIC.AvZfwd.resize(fPan.nF, bPan.nP);
    mgAIC.AwZfwd.resize(fPan.nF, bPan.nP);
    mgAIC.BuZfwd.resize(fPan.nF, bPan.nP);
    mgAIC.BvZfwd.resize(fPan.nF, bPan.nP);
    mgAIC.BwZfwd.resize(fPan.nF, bPan.nP);
    mgAIC.CuZfwd.resize(fPan.nF, fPan.nF);
    mgAIC.CvZfwd.resize(fPan.nF, fPan.nF);
    mgAIC.CwZfwd.resize(fPan.nF, fPan.nF);

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
            if (j == sp.sI[jj] && i == sp.fI[jj][ii]) {
                /// Center
                dist(0) = fPan.CG(i,0) - bPan.CG(j,0);
                dist(1) = fPan.CG(i,1) - bPan.CG(j,1);
                dist(2) = fPan.CG(i,2) - bPan.CG(j,2);
                dist = glob2loc * dist;
                trsfColloc = fPan.CG.row(i).transpose();
                trsfColloc = glob2loc * trsfColloc;

                coeffSP = split_panel(trsfCorner1(0), trsfCorner2(0), trsfCorner3(0), trsfCorner4(0),
                                      trsfCorner1(1), trsfCorner2(1), trsfCorner3(1), trsfCorner4(1),
                                      trsfColloc(0), trsfColloc(1), dist(2),
                                      bPan.l(j,0), bPan.l(j,1), bPan.l(j,2), bPan.p(j,0), bPan.p(j,1), bPan.p(j,2),
                                      bPan.n(j,0), bPan.n(j,1), bPan.n(j,2));
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

                /// X-fwd (minigrid)
                dist(0) = fPan.CG(i,0)-fPan.deltaMG - bPan.CG(j,0);
                dist(1) = fPan.CG(i,1) - bPan.CG(j,1);
                dist(2) = fPan.CG(i,2) - bPan.CG(j,2);
                dist = glob2loc * dist;
                trsfColloc = fPan.CG.row(i).transpose();
                trsfColloc(0) -= fPan.deltaMG;
                trsfColloc = glob2loc * trsfColloc;

                coeffSP = split_panel(trsfCorner1(0), trsfCorner2(0), trsfCorner3(0), trsfCorner4(0),
                                      trsfCorner1(1), trsfCorner2(1), trsfCorner3(1), trsfCorner4(1),
                                      trsfColloc(0), trsfColloc(1), dist(2),
                                      bPan.l(j,0), bPan.l(j,1), bPan.l(j,2), bPan.p(j,0), bPan.p(j,1), bPan.p(j,2),
                                      bPan.n(j,0), bPan.n(j,1), bPan.n(j,2));
                spAIC.AuXbwd[jj].row(ii) = coeffSP[0];
                spAIC.AvXbwd[jj].row(ii) = coeffSP[1];
                spAIC.AwXbwd[jj].row(ii) = coeffSP[2];
                spAIC.BuXbwd[jj].row(ii) = coeffSP[3];
                spAIC.BvXbwd[jj].row(ii) = coeffSP[4];
                spAIC.BwXbwd[jj].row(ii) = coeffSP[5];
                // Set global AIC to 0
                mgAIC.AuXbwd(i,j) = 0;
                mgAIC.AvXbwd(i,j) = 0;
                mgAIC.AwXbwd(i,j) = 0;
                mgAIC.BuXbwd(i,j) = 0;
                mgAIC.BvXbwd(i,j) = 0;
                mgAIC.BwXbwd(i,j) = 0;

                /// X-fwd (minigrid)
                dist(0) = fPan.CG(i,0)+fPan.deltaMG - bPan.CG(j,0);
                dist(1) = fPan.CG(i,1) - bPan.CG(j,1);
                dist(2) = fPan.CG(i,2) - bPan.CG(j,2);
                dist = glob2loc * dist;
                trsfColloc = fPan.CG.row(i).transpose();
                trsfColloc(0) += fPan.deltaMG;
                trsfColloc = glob2loc * trsfColloc;

                coeffSP = split_panel(trsfCorner1(0), trsfCorner2(0), trsfCorner3(0), trsfCorner4(0),
                                      trsfCorner1(1), trsfCorner2(1), trsfCorner3(1), trsfCorner4(1),
                                      trsfColloc(0), trsfColloc(1), dist(2),
                                      bPan.l(j,0), bPan.l(j,1), bPan.l(j,2), bPan.p(j,0), bPan.p(j,1), bPan.p(j,2),
                                      bPan.n(j,0), bPan.n(j,1), bPan.n(j,2));
                spAIC.AuXfwd[jj].row(ii) = coeffSP[0];
                spAIC.AvXfwd[jj].row(ii) = coeffSP[1];
                spAIC.AwXfwd[jj].row(ii) = coeffSP[2];
                spAIC.BuXfwd[jj].row(ii) = coeffSP[3];
                spAIC.BvXfwd[jj].row(ii) = coeffSP[4];
                spAIC.BwXfwd[jj].row(ii) = coeffSP[5];
                // Set global AIC to 0
                mgAIC.AuXfwd(i,j) = 0;
                mgAIC.AvXfwd(i,j) = 0;
                mgAIC.AwXfwd(i,j) = 0;
                mgAIC.BuXfwd(i,j) = 0;
                mgAIC.BvXfwd(i,j) = 0;
                mgAIC.BwXfwd(i,j) = 0;

                /// Y-bwd (minigrid)
                dist(0) = fPan.CG(i,0) - bPan.CG(j,0);
                dist(1) = fPan.CG(i,1)-fPan.deltaMG - bPan.CG(j,1);
                dist(2) = fPan.CG(i,2) - bPan.CG(j,2);
                dist = glob2loc * dist;
                trsfColloc = fPan.CG.row(i).transpose();
                trsfColloc(1) -= fPan.deltaMG;
                trsfColloc = glob2loc * trsfColloc;

                coeffSP = split_panel(trsfCorner1(0), trsfCorner2(0), trsfCorner3(0), trsfCorner4(0),
                                      trsfCorner1(1), trsfCorner2(1), trsfCorner3(1), trsfCorner4(1),
                                      trsfColloc(0), trsfColloc(1), dist(2),
                                      bPan.l(j,0), bPan.l(j,1), bPan.l(j,2), bPan.p(j,0), bPan.p(j,1), bPan.p(j,2),
                                      bPan.n(j,0), bPan.n(j,1), bPan.n(j,2));
                spAIC.AuYbwd[jj].row(ii) = coeffSP[0];
                spAIC.AvYbwd[jj].row(ii) = coeffSP[1];
                spAIC.AwYbwd[jj].row(ii) = coeffSP[2];
                spAIC.BuYbwd[jj].row(ii) = coeffSP[3];
                spAIC.BvYbwd[jj].row(ii) = coeffSP[4];
                spAIC.BwYbwd[jj].row(ii) = coeffSP[5];
                // Set global AIC to 0
                mgAIC.AuYbwd(i,j) = 0;
                mgAIC.AvYbwd(i,j) = 0;
                mgAIC.AwYbwd(i,j) = 0;
                mgAIC.BuYbwd(i,j) = 0;
                mgAIC.BvYbwd(i,j) = 0;
                mgAIC.BwYbwd(i,j) = 0;

                /// Y-fwd (minigrid)
                dist(0) = fPan.CG(i,0) - bPan.CG(j,0);
                dist(1) = fPan.CG(i,1)+fPan.deltaMG - bPan.CG(j,1);
                dist(2) = fPan.CG(i,2) - bPan.CG(j,2);
                dist = glob2loc * dist;
                trsfColloc = fPan.CG.row(i).transpose();
                trsfColloc(1) += fPan.deltaMG;
                trsfColloc = glob2loc * trsfColloc;

                coeffSP = split_panel(trsfCorner1(0), trsfCorner2(0), trsfCorner3(0), trsfCorner4(0),
                                      trsfCorner1(1), trsfCorner2(1), trsfCorner3(1), trsfCorner4(1),
                                      trsfColloc(0), trsfColloc(1), dist(2),
                                      bPan.l(j,0), bPan.l(j,1), bPan.l(j,2), bPan.p(j,0), bPan.p(j,1), bPan.p(j,2),
                                      bPan.n(j,0), bPan.n(j,1), bPan.n(j,2));
                spAIC.AuYfwd[jj].row(ii) = coeffSP[0];
                spAIC.AvYfwd[jj].row(ii) = coeffSP[1];
                spAIC.AwYfwd[jj].row(ii) = coeffSP[2];
                spAIC.BuYfwd[jj].row(ii) = coeffSP[3];
                spAIC.BvYfwd[jj].row(ii) = coeffSP[4];
                spAIC.BwYfwd[jj].row(ii) = coeffSP[5];
                // Set global AIC to 0
                mgAIC.AuYfwd(i,j) = 0;
                mgAIC.AvYfwd(i,j) = 0;
                mgAIC.AwYfwd(i,j) = 0;
                mgAIC.BuYfwd(i,j) = 0;
                mgAIC.BvYfwd(i,j) = 0;
                mgAIC.BwYfwd(i,j) = 0;

                /// Z-bwd (minigrid)
                dist(0) = fPan.CG(i,0) - bPan.CG(j,0);
                dist(1) = fPan.CG(i,1) - bPan.CG(j,1);
                dist(2) = fPan.CG(i,2)-fPan.deltaMG - bPan.CG(j,2);
                dist = glob2loc * dist;
                trsfColloc = fPan.CG.row(i).transpose();
                trsfColloc(2) -= fPan.deltaMG;
                trsfColloc = glob2loc * trsfColloc;

                coeffSP = split_panel(trsfCorner1(0), trsfCorner2(0), trsfCorner3(0), trsfCorner4(0),
                                      trsfCorner1(1), trsfCorner2(1), trsfCorner3(1), trsfCorner4(1),
                                      trsfColloc(0), trsfColloc(1), dist(2),
                                      bPan.l(j,0), bPan.l(j,1), bPan.l(j,2), bPan.p(j,0), bPan.p(j,1), bPan.p(j,2),
                                      bPan.n(j,0), bPan.n(j,1), bPan.n(j,2));
                spAIC.AuZbwd[jj].row(ii) = coeffSP[0];
                spAIC.AvZbwd[jj].row(ii) = coeffSP[1];
                spAIC.AwZbwd[jj].row(ii) = coeffSP[2];
                spAIC.BuZbwd[jj].row(ii) = coeffSP[3];
                spAIC.BvZbwd[jj].row(ii) = coeffSP[4];
                spAIC.BwZbwd[jj].row(ii) = coeffSP[5];
                // Set global AIC to 0
                mgAIC.AuZbwd(i,j) = 0;
                mgAIC.AvZbwd(i,j) = 0;
                mgAIC.AwZbwd(i,j) = 0;
                mgAIC.BuZbwd(i,j) = 0;
                mgAIC.BvZbwd(i,j) = 0;
                mgAIC.BwZbwd(i,j) = 0;

                /// Z-fwd (minigrid)
                dist(0) = fPan.CG(i,0) - bPan.CG(j,0);
                dist(1) = fPan.CG(i,1) - bPan.CG(j,1);
                dist(2) = fPan.CG(i,2)+fPan.deltaMG - bPan.CG(j,2);
                dist = glob2loc * dist;
                trsfColloc = fPan.CG.row(i).transpose();
                trsfColloc(2) += fPan.deltaMG;
                trsfColloc = glob2loc * trsfColloc;

                coeffSP = split_panel(trsfCorner1(0), trsfCorner2(0), trsfCorner3(0), trsfCorner4(0),
                                      trsfCorner1(1), trsfCorner2(1), trsfCorner3(1), trsfCorner4(1),
                                      trsfColloc(0), trsfColloc(1), dist(2),
                                      bPan.l(j,0), bPan.l(j,1), bPan.l(j,2), bPan.p(j,0), bPan.p(j,1), bPan.p(j,2),
                                      bPan.n(j,0), bPan.n(j,1), bPan.n(j,2));
                spAIC.AuZfwd[jj].row(ii) = coeffSP[0];
                spAIC.AvZfwd[jj].row(ii) = coeffSP[1];
                spAIC.AwZfwd[jj].row(ii) = coeffSP[2];
                spAIC.BuZfwd[jj].row(ii) = coeffSP[3];
                spAIC.BvZfwd[jj].row(ii) = coeffSP[4];
                spAIC.BwZfwd[jj].row(ii) = coeffSP[5];
                // Set global AIC to 0
                mgAIC.AuZfwd(i,j) = 0;
                mgAIC.AvZfwd(i,j) = 0;
                mgAIC.AwZfwd(i,j) = 0;
                mgAIC.BuZfwd(i,j) = 0;
                mgAIC.BvZfwd(i,j) = 0;
                mgAIC.BwZfwd(i,j) = 0;

                // Set flags
                if (ii == 0)
                    FLAG = 1;
                if (ii < sp.fI[jj].size())
                    i++;
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

                /// X-bwd (minigrid)
                dist(0) = fPan.CG(i,0)-fPan.deltaMG - bPan.CG(j,0);
                dist(1) = fPan.CG(i,1) - bPan.CG(j,1);
                dist(2) = fPan.CG(i,2) - bPan.CG(j,2);
                dist = glob2loc * dist;
                trsfColloc = fPan.CG.row(i).transpose();
                trsfColloc(0) -= fPan.deltaMG;
                trsfColloc = glob2loc * trsfColloc;

                coeffBF = infcBF(0, trsfCorner1(0), trsfCorner2(0), trsfCorner3(0), trsfCorner4(0),
                                 trsfCorner1(1), trsfCorner2(1), trsfCorner3(1), trsfCorner4(1),
                                 trsfColloc(0), trsfColloc(1), dist(2),
                                 bPan.l(j, 0), bPan.l(j, 1), bPan.l(j, 2), bPan.p(j, 0), bPan.p(j, 1), bPan.p(j, 2),
                                 bPan.n(j, 0), bPan.n(j, 1), bPan.n(j, 2));
                mgAIC.AuXbwd(i, j) = coeffBF[0];
                mgAIC.AvXbwd(i, j) = coeffBF[1];
                mgAIC.AwXbwd(i, j) = coeffBF[2];
                mgAIC.BuXbwd(i, j) = coeffBF[3];
                mgAIC.BvXbwd(i, j) = coeffBF[4];
                mgAIC.BwXbwd(i, j) = coeffBF[5];

                // X-fwd (minigrid)
                dist(0) = fPan.CG(i,0)+fPan.deltaMG - bPan.CG(j,0);
                dist(1) = fPan.CG(i,1) - bPan.CG(j,1);
                dist(2) = fPan.CG(i,2) - bPan.CG(j,2);
                dist = glob2loc * dist;
                trsfColloc = fPan.CG.row(i).transpose();
                trsfColloc(0) += fPan.deltaMG;
                trsfColloc = glob2loc * trsfColloc;

                coeffBF = infcBF(0, trsfCorner1(0) ,trsfCorner2(0), trsfCorner3(0), trsfCorner4(0),
                                 trsfCorner1(1), trsfCorner2(1), trsfCorner3(1), trsfCorner4(1),
                                 trsfColloc(0), trsfColloc(1), dist(2),
                                 bPan.l(j,0), bPan.l(j,1), bPan.l(j,2), bPan.p(j,0), bPan.p(j,1), bPan.p(j,2),
                                 bPan.n(j,0), bPan.n(j,1), bPan.n(j,2));
                mgAIC.AuXfwd(i,j) = coeffBF[0];
                mgAIC.AvXfwd(i,j) = coeffBF[1];
                mgAIC.AwXfwd(i,j) = coeffBF[2];
                mgAIC.BuXfwd(i,j) = coeffBF[3];
                mgAIC.BvXfwd(i,j) = coeffBF[4];
                mgAIC.BwXfwd(i,j) = coeffBF[5];

                // Y-bwd (minigrid)
                dist(0) = fPan.CG(i,0) - bPan.CG(j,0);
                dist(1) = fPan.CG(i,1)-fPan.deltaMG - bPan.CG(j,1);
                dist(2) = fPan.CG(i,2) - bPan.CG(j,2);
                dist = glob2loc * dist;
                trsfColloc = fPan.CG.row(i).transpose();
                trsfColloc(1) -= fPan.deltaMG;
                trsfColloc = glob2loc * trsfColloc;

                coeffBF = infcBF(0, trsfCorner1(0) ,trsfCorner2(0), trsfCorner3(0), trsfCorner4(0),
                                 trsfCorner1(1), trsfCorner2(1), trsfCorner3(1), trsfCorner4(1),
                                 trsfColloc(0), trsfColloc(1), dist(2),
                                 bPan.l(j,0), bPan.l(j,1), bPan.l(j,2), bPan.p(j,0), bPan.p(j,1), bPan.p(j,2),
                                 bPan.n(j,0), bPan.n(j,1), bPan.n(j,2));
                mgAIC.AuYbwd(i,j) = coeffBF[0];
                mgAIC.AvYbwd(i,j) = coeffBF[1];
                mgAIC.AwYbwd(i,j) = coeffBF[2];
                mgAIC.BuYbwd(i,j) = coeffBF[3];
                mgAIC.BvYbwd(i,j) = coeffBF[4];
                mgAIC.BwYbwd(i,j) = coeffBF[5];

                // Y-fwd (minigrid)
                dist(0) = fPan.CG(i,0) - bPan.CG(j,0);
                dist(1) = fPan.CG(i,1)+fPan.deltaMG - bPan.CG(j,1);
                dist(2) = fPan.CG(i,2) - bPan.CG(j,2);
                dist = glob2loc * dist;
                trsfColloc = fPan.CG.row(i).transpose();
                trsfColloc(1) += fPan.deltaMG;
                trsfColloc = glob2loc * trsfColloc;

                coeffBF = infcBF(0, trsfCorner1(0) ,trsfCorner2(0), trsfCorner3(0), trsfCorner4(0),
                                 trsfCorner1(1), trsfCorner2(1), trsfCorner3(1), trsfCorner4(1),
                                 trsfColloc(0), trsfColloc(1), dist(2),
                                 bPan.l(j,0), bPan.l(j,1), bPan.l(j,2), bPan.p(j,0), bPan.p(j,1), bPan.p(j,2),
                                 bPan.n(j,0), bPan.n(j,1), bPan.n(j,2));
                mgAIC.AuYfwd(i,j) = coeffBF[0];
                mgAIC.AvYfwd(i,j) = coeffBF[1];
                mgAIC.AwYfwd(i,j) = coeffBF[2];
                mgAIC.BuYfwd(i,j) = coeffBF[3];
                mgAIC.BvYfwd(i,j) = coeffBF[4];
                mgAIC.BwYfwd(i,j) = coeffBF[5];

                // Z-bwd (minigrid)
                dist(0) = fPan.CG(i,0) - bPan.CG(j,0);
                dist(1) = fPan.CG(i,1) - bPan.CG(j,1);
                dist(2) = fPan.CG(i,2)-fPan.deltaMG - bPan.CG(j,2);
                dist = glob2loc * dist;
                trsfColloc = fPan.CG.row(i).transpose();
                trsfColloc(2) -= fPan.deltaMG;
                trsfColloc = glob2loc * trsfColloc;

                coeffBF = infcBF(0, trsfCorner1(0) ,trsfCorner2(0), trsfCorner3(0), trsfCorner4(0),
                                 trsfCorner1(1), trsfCorner2(1), trsfCorner3(1), trsfCorner4(1),
                                 trsfColloc(0), trsfColloc(1), dist(2),
                                 bPan.l(j,0), bPan.l(j,1), bPan.l(j,2), bPan.p(j,0), bPan.p(j,1), bPan.p(j,2),
                                 bPan.n(j,0), bPan.n(j,1), bPan.n(j,2));
                mgAIC.AuZbwd(i,j) = coeffBF[0];
                mgAIC.AvZbwd(i,j) = coeffBF[1];
                mgAIC.AwZbwd(i,j) = coeffBF[2];
                mgAIC.BuZbwd(i,j) = coeffBF[3];
                mgAIC.BvZbwd(i,j) = coeffBF[4];
                mgAIC.BwZbwd(i,j) = coeffBF[5];

                // Z-fwd (minigrid)
                dist(0) = fPan.CG(i,0) - bPan.CG(j,0);
                dist(1) = fPan.CG(i,1) - bPan.CG(j,1);
                dist(2) = fPan.CG(i,2)+fPan.deltaMG - bPan.CG(j,2);
                dist = glob2loc * dist;
                trsfColloc = fPan.CG.row(i).transpose();
                trsfColloc(2) += fPan.deltaMG;
                trsfColloc = glob2loc * trsfColloc;

                coeffBF = infcBF(0, trsfCorner1(0) ,trsfCorner2(0), trsfCorner3(0), trsfCorner4(0),
                                 trsfCorner1(1), trsfCorner2(1), trsfCorner3(1), trsfCorner4(1),
                                 trsfColloc(0), trsfColloc(1), dist(2),
                                 bPan.l(j,0), bPan.l(j,1), bPan.l(j,2), bPan.p(j,0), bPan.p(j,1), bPan.p(j,2),
                                 bPan.n(j,0), bPan.n(j,1), bPan.n(j,2));
                mgAIC.AuZfwd(i,j) = coeffBF[0];
                mgAIC.AvZfwd(i,j) = coeffBF[1];
                mgAIC.AwZfwd(i,j) = coeffBF[2];
                mgAIC.BuZfwd(i,j) = coeffBF[3];
                mgAIC.BvZfwd(i,j) = coeffBF[4];
                mgAIC.BwZfwd(i,j) = coeffBF[5];
            }
            if (i == fPan.nF-1 && FLAG && jj < sp.sI.size())
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

            // panel to field for minigrid
            for (int i = 0; i < fPan.nF; ++i) {
                dist(0) = fPan.CG(i, 0)-fPan.deltaMG - bPan.CG(j, 0);
                dist(1) = fPan.CG(i, 1) + bPan.CG(j, 1);
                dist(2) = fPan.CG(i, 2) - bPan.CG(j, 2);
                dist = glob2loc * dist;
                trsfColloc = fPan.CG.row(i).transpose();
                trsfColloc(0) -= fPan.deltaMG;
                trsfColloc = glob2loc * trsfColloc;

                coeffBF = infcBF(0, trsfCorner1(0) ,trsfCorner2(0), trsfCorner3(0), trsfCorner4(0),
                                 trsfCorner1(1), trsfCorner2(1), trsfCorner3(1), trsfCorner4(1),
                                 trsfColloc(0), trsfColloc(1), dist(2),
                                 bPan.l(j,0), -bPan.l(j,1), bPan.l(j,2), -bPan.p(j,0), bPan.p(j,1), -bPan.p(j,2),
                                 bPan.n(j,0), -bPan.n(j,1), bPan.n(j,2));

                mgAIC.AuXbwd(i,j) -= coeffBF[0];
                mgAIC.AvXbwd(i,j) -= coeffBF[1];
                mgAIC.AwXbwd(i,j) -= coeffBF[2];
                mgAIC.BuXbwd(i,j) -= coeffBF[3];
                mgAIC.BvXbwd(i,j) -= coeffBF[4];
                mgAIC.BwXbwd(i,j) -= coeffBF[5];
            }
            for (int i = 0; i < fPan.nF; ++i) {
                dist(0) = fPan.CG(i, 0)+fPan.deltaMG - bPan.CG(j, 0);
                dist(1) = fPan.CG(i, 1) + bPan.CG(j, 1);
                dist(2) = fPan.CG(i, 2) - bPan.CG(j, 2);
                dist = glob2loc * dist;
                trsfColloc = fPan.CG.row(i).transpose();
                trsfColloc(0) += fPan.deltaMG;
                trsfColloc = glob2loc * trsfColloc;

                coeffBF = infcBF(0, trsfCorner1(0) ,trsfCorner2(0), trsfCorner3(0), trsfCorner4(0),
                                 trsfCorner1(1), trsfCorner2(1), trsfCorner3(1), trsfCorner4(1),
                                 trsfColloc(0), trsfColloc(1), dist(2),
                                 bPan.l(j,0), -bPan.l(j,1), bPan.l(j,2), -bPan.p(j,0), bPan.p(j,1), -bPan.p(j,2),
                                 bPan.n(j,0), -bPan.n(j,1), bPan.n(j,2));

                mgAIC.AuXfwd(i,j) -= coeffBF[0];
                mgAIC.AvXfwd(i,j) -= coeffBF[1];
                mgAIC.AwXfwd(i,j) -= coeffBF[2];
                mgAIC.BuXfwd(i,j) -= coeffBF[3];
                mgAIC.BvXfwd(i,j) -= coeffBF[4];
                mgAIC.BwXfwd(i,j) -= coeffBF[5];
            }
            for (int i = 0; i < fPan.nF; ++i) {
                dist(0) = fPan.CG(i, 0) - bPan.CG(j, 0);
                dist(1) = fPan.CG(i, 1)-fPan.deltaMG + bPan.CG(j, 1);
                dist(2) = fPan.CG(i, 2) - bPan.CG(j, 2);
                dist = glob2loc * dist;
                trsfColloc = fPan.CG.row(i).transpose();
                trsfColloc(1) -= fPan.deltaMG;
                trsfColloc = glob2loc * trsfColloc;

                coeffBF = infcBF(0, trsfCorner1(0) ,trsfCorner2(0), trsfCorner3(0), trsfCorner4(0),
                                 trsfCorner1(1), trsfCorner2(1), trsfCorner3(1), trsfCorner4(1),
                                 trsfColloc(0), trsfColloc(1), dist(2),
                                 bPan.l(j,0), -bPan.l(j,1), bPan.l(j,2), -bPan.p(j,0), bPan.p(j,1), -bPan.p(j,2),
                                 bPan.n(j,0), -bPan.n(j,1), bPan.n(j,2));

                mgAIC.AuYbwd(i,j) -= coeffBF[0];
                mgAIC.AvYbwd(i,j) -= coeffBF[1];
                mgAIC.AwYbwd(i,j) -= coeffBF[2];
                mgAIC.BuYbwd(i,j) -= coeffBF[3];
                mgAIC.BvYbwd(i,j) -= coeffBF[4];
                mgAIC.BwYbwd(i,j) -= coeffBF[5];
            }
            for (int i = 0; i < fPan.nF; ++i) {
                dist(0) = fPan.CG(i, 0) - bPan.CG(j, 0);
                dist(1) = fPan.CG(i, 1)+fPan.deltaMG + bPan.CG(j, 1);
                dist(2) = fPan.CG(i, 2) - bPan.CG(j, 2);
                dist = glob2loc * dist;
                trsfColloc = fPan.CG.row(i).transpose();
                trsfColloc(1) += fPan.deltaMG;
                trsfColloc = glob2loc * trsfColloc;

                coeffBF = infcBF(0, trsfCorner1(0) ,trsfCorner2(0), trsfCorner3(0), trsfCorner4(0),
                                 trsfCorner1(1), trsfCorner2(1), trsfCorner3(1), trsfCorner4(1),
                                 trsfColloc(0), trsfColloc(1), dist(2),
                                 bPan.l(j,0), -bPan.l(j,1), bPan.l(j,2), -bPan.p(j,0), bPan.p(j,1), -bPan.p(j,2),
                                 bPan.n(j,0), -bPan.n(j,1), bPan.n(j,2));

                mgAIC.AuYfwd(i,j) -= coeffBF[0];
                mgAIC.AvYfwd(i,j) -= coeffBF[1];
                mgAIC.AwYfwd(i,j) -= coeffBF[2];
                mgAIC.BuYfwd(i,j) -= coeffBF[3];
                mgAIC.BvYfwd(i,j) -= coeffBF[4];
                mgAIC.BwYfwd(i,j) -= coeffBF[5];
            }
            for (int i = 0; i < fPan.nF; ++i) {
                dist(0) = fPan.CG(i, 0) - bPan.CG(j, 0);
                dist(1) = fPan.CG(i, 1) + bPan.CG(j, 1);
                dist(2) = fPan.CG(i, 2)-fPan.deltaMG - bPan.CG(j, 2);
                dist = glob2loc * dist;
                trsfColloc = fPan.CG.row(i).transpose();
                trsfColloc(2) -= fPan.deltaMG;
                trsfColloc = glob2loc * trsfColloc;

                coeffBF = infcBF(0, trsfCorner1(0) ,trsfCorner2(0), trsfCorner3(0), trsfCorner4(0),
                                 trsfCorner1(1), trsfCorner2(1), trsfCorner3(1), trsfCorner4(1),
                                 trsfColloc(0), trsfColloc(1), dist(2),
                                 bPan.l(j,0), -bPan.l(j,1), bPan.l(j,2), -bPan.p(j,0), bPan.p(j,1), -bPan.p(j,2),
                                 bPan.n(j,0), -bPan.n(j,1), bPan.n(j,2));

                mgAIC.AuZbwd(i,j) -= coeffBF[0];
                mgAIC.AvZbwd(i,j) -= coeffBF[1];
                mgAIC.AwZbwd(i,j) -= coeffBF[2];
                mgAIC.BuZbwd(i,j) -= coeffBF[3];
                mgAIC.BvZbwd(i,j) -= coeffBF[4];
                mgAIC.BwZbwd(i,j) -= coeffBF[5];
            }
            for (int i = 0; i < fPan.nF; ++i) {
                dist(0) = fPan.CG(i, 0) - bPan.CG(j, 0);
                dist(1) = fPan.CG(i, 1) + bPan.CG(j, 1);
                dist(2) = fPan.CG(i, 2)+fPan.deltaMG - bPan.CG(j, 2);
                dist = glob2loc * dist;
                trsfColloc = fPan.CG.row(i).transpose();
                trsfColloc(2) += fPan.deltaMG;
                trsfColloc = glob2loc * trsfColloc;

                coeffBF = infcBF(0, trsfCorner1(0) ,trsfCorner2(0), trsfCorner3(0), trsfCorner4(0),
                                 trsfCorner1(1), trsfCorner2(1), trsfCorner3(1), trsfCorner4(1),
                                 trsfColloc(0), trsfColloc(1), dist(2),
                                 bPan.l(j,0), -bPan.l(j,1), bPan.l(j,2), -bPan.p(j,0), bPan.p(j,1), -bPan.p(j,2),
                                 bPan.n(j,0), -bPan.n(j,1), bPan.n(j,2));

                mgAIC.AuZfwd(i,j) -= coeffBF[0];
                mgAIC.AvZfwd(i,j) -= coeffBF[1];
                mgAIC.AwZfwd(i,j) -= coeffBF[2];
                mgAIC.BuZfwd(i,j) -= coeffBF[3];
                mgAIC.BvZfwd(i,j) -= coeffBF[4];
                mgAIC.BwZfwd(i,j) -= coeffBF[5];
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

        // panel to field for minigrid
        for (int i = 0; i < fPan.nF; ++i) {
            dist(0) = fPan.CG(i, 0)-fPan.deltaMG - wPan.CG(j, 0);
            dist(1) = fPan.CG(i, 1) - wPan.CG(j, 1);
            dist(2) = fPan.CG(i, 2) - wPan.CG(j, 2);
            dist = glob2loc * dist;
            trsfColloc = fPan.CG.row(i).transpose();
            trsfColloc(0) -= fPan.deltaMG;
            trsfColloc = glob2loc * trsfColloc;

            coeffBF = infcBF(1, trsfCorner1(0) ,trsfCorner2(0), trsfCorner3(0), trsfCorner4(0),
                             trsfCorner1(1), trsfCorner2(1), trsfCorner3(1), trsfCorner4(1),
                             trsfColloc(0), trsfColloc(1), dist(2),
                             wPan.l(j,0), wPan.l(j,1), wPan.l(j,2), wPan.p(j,0), wPan.p(j,1), wPan.p(j,2),
                             wPan.n(j,0), wPan.n(j,1), wPan.n(j,2));

            mgAIC.AuXbwd(i, j * bPan.nC_) += coeffBF[0];
            mgAIC.AuXbwd(i, (j + 1) * bPan.nC_ - 1) -= coeffBF[0];
            mgAIC.AvXbwd(i, j * bPan.nC_) += coeffBF[1];
            mgAIC.AvXbwd(i, (j + 1) * bPan.nC_ - 1) -= coeffBF[1];
            mgAIC.AwXbwd(i, j * bPan.nC_) += coeffBF[2];
            mgAIC.AwXbwd(i, (j + 1) * bPan.nC_ - 1) -= coeffBF[2];
        }
        for (int i = 0; i < fPan.nF; ++i) {
            dist(0) = fPan.CG(i, 0)+fPan.deltaMG - wPan.CG(j, 0);
            dist(1) = fPan.CG(i, 1) - wPan.CG(j, 1);
            dist(2) = fPan.CG(i, 2) - wPan.CG(j, 2);
            dist = glob2loc * dist;
            trsfColloc = fPan.CG.row(i).transpose();
            trsfColloc(0) += fPan.deltaMG;
            trsfColloc = glob2loc * trsfColloc;

            coeffBF = infcBF(1, trsfCorner1(0) ,trsfCorner2(0), trsfCorner3(0), trsfCorner4(0),
                             trsfCorner1(1), trsfCorner2(1), trsfCorner3(1), trsfCorner4(1),
                             trsfColloc(0), trsfColloc(1), dist(2),
                             wPan.l(j,0), wPan.l(j,1), wPan.l(j,2), wPan.p(j,0), wPan.p(j,1), wPan.p(j,2),
                             wPan.n(j,0), wPan.n(j,1), wPan.n(j,2));

            mgAIC.AuXfwd(i, j * bPan.nC_) += coeffBF[0];
            mgAIC.AuXfwd(i, (j + 1) * bPan.nC_ - 1) -= coeffBF[0];
            mgAIC.AvXfwd(i, j * bPan.nC_) += coeffBF[1];
            mgAIC.AvXfwd(i, (j + 1) * bPan.nC_ - 1) -= coeffBF[1];
            mgAIC.AwXfwd(i, j * bPan.nC_) += coeffBF[2];
            mgAIC.AwXfwd(i, (j + 1) * bPan.nC_ - 1) -= coeffBF[2];
        }
        for (int i = 0; i < fPan.nF; ++i) {
            dist(0) = fPan.CG(i, 0) - wPan.CG(j, 0);
            dist(1) = fPan.CG(i, 1)-fPan.deltaMG - wPan.CG(j, 1);
            dist(2) = fPan.CG(i, 2) - wPan.CG(j, 2);
            dist = glob2loc * dist;
            trsfColloc = fPan.CG.row(i).transpose();
            trsfColloc(1) -= fPan.deltaMG;
            trsfColloc = glob2loc * trsfColloc;

            coeffBF = infcBF(1, trsfCorner1(0) ,trsfCorner2(0), trsfCorner3(0), trsfCorner4(0),
                             trsfCorner1(1), trsfCorner2(1), trsfCorner3(1), trsfCorner4(1),
                             trsfColloc(0), trsfColloc(1), dist(2),
                             wPan.l(j,0), wPan.l(j,1), wPan.l(j,2), wPan.p(j,0), wPan.p(j,1), wPan.p(j,2),
                             wPan.n(j,0), wPan.n(j,1), wPan.n(j,2));

            mgAIC.AuYbwd(i, j * bPan.nC_) += coeffBF[0];
            mgAIC.AuYbwd(i, (j + 1) * bPan.nC_ - 1) -= coeffBF[0];
            mgAIC.AvYbwd(i, j * bPan.nC_) += coeffBF[1];
            mgAIC.AvYbwd(i, (j + 1) * bPan.nC_ - 1) -= coeffBF[1];
            mgAIC.AwYbwd(i, j * bPan.nC_) += coeffBF[2];
            mgAIC.AwYbwd(i, (j + 1) * bPan.nC_ - 1) -= coeffBF[2];
        }
        for (int i = 0; i < fPan.nF; ++i) {
            dist(0) = fPan.CG(i, 0) - wPan.CG(j, 0);
            dist(1) = fPan.CG(i, 1)+fPan.deltaMG - wPan.CG(j, 1);
            dist(2) = fPan.CG(i, 2) - wPan.CG(j, 2);
            dist = glob2loc * dist;
            trsfColloc = fPan.CG.row(i).transpose();
            trsfColloc(1) += fPan.deltaMG;
            trsfColloc = glob2loc * trsfColloc;

            coeffBF = infcBF(1, trsfCorner1(0) ,trsfCorner2(0), trsfCorner3(0), trsfCorner4(0),
                             trsfCorner1(1), trsfCorner2(1), trsfCorner3(1), trsfCorner4(1),
                             trsfColloc(0), trsfColloc(1), dist(2),
                             wPan.l(j,0), wPan.l(j,1), wPan.l(j,2), wPan.p(j,0), wPan.p(j,1), wPan.p(j,2),
                             wPan.n(j,0), wPan.n(j,1), wPan.n(j,2));

            mgAIC.AuYfwd(i, j * bPan.nC_) += coeffBF[0];
            mgAIC.AuYfwd(i, (j + 1) * bPan.nC_ - 1) -= coeffBF[0];
            mgAIC.AvYfwd(i, j * bPan.nC_) += coeffBF[1];
            mgAIC.AvYfwd(i, (j + 1) * bPan.nC_ - 1) -= coeffBF[1];
            mgAIC.AwYfwd(i, j * bPan.nC_) += coeffBF[2];
            mgAIC.AwYfwd(i, (j + 1) * bPan.nC_ - 1) -= coeffBF[2];
        }
        for (int i = 0; i < fPan.nF; ++i) {
            dist(0) = fPan.CG(i, 0) - wPan.CG(j, 0);
            dist(1) = fPan.CG(i, 1) - wPan.CG(j, 1);
            dist(2) = fPan.CG(i, 2)-fPan.deltaMG - wPan.CG(j, 2);
            dist = glob2loc * dist;
            trsfColloc = fPan.CG.row(i).transpose();
            trsfColloc(2) -= fPan.deltaMG;
            trsfColloc = glob2loc * trsfColloc;

            coeffBF = infcBF(1, trsfCorner1(0) ,trsfCorner2(0), trsfCorner3(0), trsfCorner4(0),
                             trsfCorner1(1), trsfCorner2(1), trsfCorner3(1), trsfCorner4(1),
                             trsfColloc(0), trsfColloc(1), dist(2),
                             wPan.l(j,0), wPan.l(j,1), wPan.l(j,2), wPan.p(j,0), wPan.p(j,1), wPan.p(j,2),
                             wPan.n(j,0), wPan.n(j,1), wPan.n(j,2));

            mgAIC.AuZbwd(i, j * bPan.nC_) += coeffBF[0];
            mgAIC.AuZbwd(i, (j + 1) * bPan.nC_ - 1) -= coeffBF[0];
            mgAIC.AvZbwd(i, j * bPan.nC_) += coeffBF[1];
            mgAIC.AvZbwd(i, (j + 1) * bPan.nC_ - 1) -= coeffBF[1];
            mgAIC.AwZbwd(i, j * bPan.nC_) += coeffBF[2];
            mgAIC.AwZbwd(i, (j + 1) * bPan.nC_ - 1) -= coeffBF[2];
        }
        for (int i = 0; i < fPan.nF; ++i) {
            dist(0) = fPan.CG(i, 0) - wPan.CG(j, 0);
            dist(1) = fPan.CG(i, 1) - wPan.CG(j, 1);
            dist(2) = fPan.CG(i, 2)+fPan.deltaMG - wPan.CG(j, 2);
            dist = glob2loc * dist;
            trsfColloc = fPan.CG.row(i).transpose();
            trsfColloc(2) += fPan.deltaMG;
            trsfColloc = glob2loc * trsfColloc;

            coeffBF = infcBF(1, trsfCorner1(0) ,trsfCorner2(0), trsfCorner3(0), trsfCorner4(0),
                             trsfCorner1(1), trsfCorner2(1), trsfCorner3(1), trsfCorner4(1),
                             trsfColloc(0), trsfColloc(1), dist(2),
                             wPan.l(j,0), wPan.l(j,1), wPan.l(j,2), wPan.p(j,0), wPan.p(j,1), wPan.p(j,2),
                             wPan.n(j,0), wPan.n(j,1), wPan.n(j,2));

            mgAIC.AuZfwd(i, j * bPan.nC_) += coeffBF[0];
            mgAIC.AuZfwd(i, (j + 1) * bPan.nC_ - 1) -= coeffBF[0];
            mgAIC.AvZfwd(i, j * bPan.nC_) += coeffBF[1];
            mgAIC.AvZfwd(i, (j + 1) * bPan.nC_ - 1) -= coeffBF[1];
            mgAIC.AwZfwd(i, j * bPan.nC_) += coeffBF[2];
            mgAIC.AwZfwd(i, (j + 1) * bPan.nC_ - 1) -= coeffBF[2];
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

            // panel to field for minigrid
            for (int i = 0; i < fPan.nF; ++i) {
                dist(0) = fPan.CG(i, 0)-fPan.deltaMG - wPan.CG(j, 0);
                dist(1) = fPan.CG(i, 1) + wPan.CG(j, 1);
                dist(2) = fPan.CG(i, 2) - wPan.CG(j, 2);
                dist = glob2loc * dist;
                trsfColloc = fPan.CG.row(i).transpose();
                trsfColloc(0) -= fPan.deltaMG;
                trsfColloc = glob2loc * trsfColloc;

                coeffBF = infcBF(1, trsfCorner1(0) ,trsfCorner2(0), trsfCorner3(0), trsfCorner4(0),
                                 trsfCorner1(1), trsfCorner2(1), trsfCorner3(1), trsfCorner4(1),
                                 trsfColloc(0), trsfColloc(1), dist(2),
                                 wPan.l(j,0), -wPan.l(j,1), wPan.l(j,2), -wPan.p(j,0), wPan.p(j,1), -wPan.p(j,2),
                                 wPan.n(j,0), -wPan.n(j,1), wPan.n(j,2));

                mgAIC.AuXbwd(i, j * bPan.nC_) -= coeffBF[0];
                mgAIC.AuXbwd(i, (j + 1) * bPan.nC_ - 1) += coeffBF[0];
                mgAIC.AvXbwd(i, j * bPan.nC_) -= coeffBF[1];
                mgAIC.AvXbwd(i, (j + 1) * bPan.nC_ - 1) += coeffBF[1];
                mgAIC.AwXbwd(i, j * bPan.nC_) -= coeffBF[2];
                mgAIC.AwXbwd(i, (j + 1) * bPan.nC_ - 1) += coeffBF[2];
            }
            for (int i = 0; i < fPan.nF; ++i) {
                dist(0) = fPan.CG(i, 0)+fPan.deltaMG - wPan.CG(j, 0);
                dist(1) = fPan.CG(i, 1) + wPan.CG(j, 1);
                dist(2) = fPan.CG(i, 2) - wPan.CG(j, 2);
                dist = glob2loc * dist;
                trsfColloc = fPan.CG.row(i).transpose();
                trsfColloc(0) += fPan.deltaMG;
                trsfColloc = glob2loc * trsfColloc;

                coeffBF = infcBF(1, trsfCorner1(0) ,trsfCorner2(0), trsfCorner3(0), trsfCorner4(0),
                                 trsfCorner1(1), trsfCorner2(1), trsfCorner3(1), trsfCorner4(1),
                                 trsfColloc(0), trsfColloc(1), dist(2),
                                 wPan.l(j,0), -wPan.l(j,1), wPan.l(j,2), -wPan.p(j,0), wPan.p(j,1), -wPan.p(j,2),
                                 wPan.n(j,0), -wPan.n(j,1), wPan.n(j,2));

                mgAIC.AuXfwd(i, j * bPan.nC_) -= coeffBF[0];
                mgAIC.AuXfwd(i, (j + 1) * bPan.nC_ - 1) += coeffBF[0];
                mgAIC.AvXfwd(i, j * bPan.nC_) -= coeffBF[1];
                mgAIC.AvXfwd(i, (j + 1) * bPan.nC_ - 1) += coeffBF[1];
                mgAIC.AwXfwd(i, j * bPan.nC_) -= coeffBF[2];
                mgAIC.AwXfwd(i, (j + 1) * bPan.nC_ - 1) += coeffBF[2];
            }
            for (int i = 0; i < fPan.nF; ++i) {
                dist(0) = fPan.CG(i, 0) - wPan.CG(j, 0);
                dist(1) = fPan.CG(i, 1)-fPan.deltaMG + wPan.CG(j, 1);
                dist(2) = fPan.CG(i, 2) - wPan.CG(j, 2);
                dist = glob2loc * dist;
                trsfColloc = fPan.CG.row(i).transpose();
                trsfColloc(1) -= fPan.deltaMG;
                trsfColloc = glob2loc * trsfColloc;

                coeffBF = infcBF(1, trsfCorner1(0) ,trsfCorner2(0), trsfCorner3(0), trsfCorner4(0),
                                 trsfCorner1(1), trsfCorner2(1), trsfCorner3(1), trsfCorner4(1),
                                 trsfColloc(0), trsfColloc(1), dist(2),
                                 wPan.l(j,0), -wPan.l(j,1), wPan.l(j,2), -wPan.p(j,0), wPan.p(j,1), -wPan.p(j,2),
                                 wPan.n(j,0), -wPan.n(j,1), wPan.n(j,2));

                mgAIC.AuYbwd(i, j * bPan.nC_) -= coeffBF[0];
                mgAIC.AuYbwd(i, (j + 1) * bPan.nC_ - 1) += coeffBF[0];
                mgAIC.AvYbwd(i, j * bPan.nC_) -= coeffBF[1];
                mgAIC.AvYbwd(i, (j + 1) * bPan.nC_ - 1) += coeffBF[1];
                mgAIC.AwYbwd(i, j * bPan.nC_) -= coeffBF[2];
                mgAIC.AwYbwd(i, (j + 1) * bPan.nC_ - 1) += coeffBF[2];
            }
            for (int i = 0; i < fPan.nF; ++i) {
                dist(0) = fPan.CG(i, 0) - wPan.CG(j, 0);
                dist(1) = fPan.CG(i, 1)+fPan.deltaMG + wPan.CG(j, 1);
                dist(2) = fPan.CG(i, 2) - wPan.CG(j, 2);
                dist = glob2loc * dist;
                trsfColloc = fPan.CG.row(i).transpose();
                trsfColloc(1) += fPan.deltaMG;
                trsfColloc = glob2loc * trsfColloc;

                coeffBF = infcBF(1, trsfCorner1(0) ,trsfCorner2(0), trsfCorner3(0), trsfCorner4(0),
                                 trsfCorner1(1), trsfCorner2(1), trsfCorner3(1), trsfCorner4(1),
                                 trsfColloc(0), trsfColloc(1), dist(2),
                                 wPan.l(j,0), -wPan.l(j,1), wPan.l(j,2), -wPan.p(j,0), wPan.p(j,1), -wPan.p(j,2),
                                 wPan.n(j,0), -wPan.n(j,1), wPan.n(j,2));

                mgAIC.AuYfwd(i, j * bPan.nC_) -= coeffBF[0];
                mgAIC.AuYfwd(i, (j + 1) * bPan.nC_ - 1) += coeffBF[0];
                mgAIC.AvYfwd(i, j * bPan.nC_) -= coeffBF[1];
                mgAIC.AvYfwd(i, (j + 1) * bPan.nC_ - 1) += coeffBF[1];
                mgAIC.AwYfwd(i, j * bPan.nC_) -= coeffBF[2];
                mgAIC.AwYfwd(i, (j + 1) * bPan.nC_ - 1) += coeffBF[2];
            }
            for (int i = 0; i < fPan.nF; ++i) {
                dist(0) = fPan.CG(i, 0) - wPan.CG(j, 0);
                dist(1) = fPan.CG(i, 1) + wPan.CG(j, 1);
                dist(2) = fPan.CG(i, 2)-fPan.deltaMG - wPan.CG(j, 2);
                dist = glob2loc * dist;
                trsfColloc = fPan.CG.row(i).transpose();
                trsfColloc(2) -= fPan.deltaMG;
                trsfColloc = glob2loc * trsfColloc;

                coeffBF = infcBF(1, trsfCorner1(0) ,trsfCorner2(0), trsfCorner3(0), trsfCorner4(0),
                                 trsfCorner1(1), trsfCorner2(1), trsfCorner3(1), trsfCorner4(1),
                                 trsfColloc(0), trsfColloc(1), dist(2),
                                 wPan.l(j,0), -wPan.l(j,1), wPan.l(j,2), -wPan.p(j,0), wPan.p(j,1), -wPan.p(j,2),
                                 wPan.n(j,0), -wPan.n(j,1), wPan.n(j,2));

                mgAIC.AuZbwd(i, j * bPan.nC_) -= coeffBF[0];
                mgAIC.AuZbwd(i, (j + 1) * bPan.nC_ - 1) += coeffBF[0];
                mgAIC.AvZbwd(i, j * bPan.nC_) -= coeffBF[1];
                mgAIC.AvZbwd(i, (j + 1) * bPan.nC_ - 1) += coeffBF[1];
                mgAIC.AwZbwd(i, j * bPan.nC_) -= coeffBF[2];
                mgAIC.AwZbwd(i, (j + 1) * bPan.nC_ - 1) += coeffBF[2];
            }
            for (int i = 0; i < fPan.nF; ++i) {
                dist(0) = fPan.CG(i, 0) - wPan.CG(j, 0);
                dist(1) = fPan.CG(i, 1) + wPan.CG(j, 1);
                dist(2) = fPan.CG(i, 2)+fPan.deltaMG - wPan.CG(j, 2);
                dist = glob2loc * dist;
                trsfColloc = fPan.CG.row(i).transpose();
                trsfColloc(2) += fPan.deltaMG;
                trsfColloc = glob2loc * trsfColloc;

                coeffBF = infcBF(1, trsfCorner1(0) ,trsfCorner2(0), trsfCorner3(0), trsfCorner4(0),
                                 trsfCorner1(1), trsfCorner2(1), trsfCorner3(1), trsfCorner4(1),
                                 trsfColloc(0), trsfColloc(1), dist(2),
                                 wPan.l(j,0), -wPan.l(j,1), wPan.l(j,2), -wPan.p(j,0), wPan.p(j,1), -wPan.p(j,2),
                                 wPan.n(j,0), -wPan.n(j,1), wPan.n(j,2));

                mgAIC.AuZfwd(i, j * bPan.nC_) -= coeffBF[0];
                mgAIC.AuZfwd(i, (j + 1) * bPan.nC_ - 1) += coeffBF[0];
                mgAIC.AvZfwd(i, j * bPan.nC_) -= coeffBF[1];
                mgAIC.AvZfwd(i, (j + 1) * bPan.nC_ - 1) += coeffBF[1];
                mgAIC.AwZfwd(i, j * bPan.nC_) -= coeffBF[2];
                mgAIC.AwZfwd(i, (j + 1) * bPan.nC_ - 1) += coeffBF[2];
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

        // Field to field for minigrid
        for (int i = 0; i < fPan.nF; ++i) {
            x = fPan.CG(i,0)-fPan.deltaMG - fPan.CG(j,0);
            y = fPan.CG(i,1) - fPan.CG(j,1);
            z = fPan.CG(i,2) - fPan.CG(j,2);
            coeffF = infcF(x, y, z, xC, yC, zC);
            mgAIC.CuXbwd(i,j) = coeffF[0];
            mgAIC.CvXbwd(i,j) = coeffF[1];
            mgAIC.CwXbwd(i,j) = coeffF[2];
        }
        for (int i = 0; i < fPan.nF; ++i) {
            x = fPan.CG(i,0)+fPan.deltaMG - fPan.CG(j,0);
            y = fPan.CG(i,1) - fPan.CG(j,1);
            z = fPan.CG(i,2) - fPan.CG(j,2);
            coeffF = infcF(x, y, z, xC, yC, zC);
            mgAIC.CuXfwd(i,j) = coeffF[0];
            mgAIC.CvXfwd(i,j) = coeffF[1];
            mgAIC.CwXfwd(i,j) = coeffF[2];
        }
        for (int i = 0; i < fPan.nF; ++i) {
            x = fPan.CG(i,0) - fPan.CG(j,0);
            y = fPan.CG(i,1)-fPan.deltaMG - fPan.CG(j,1);
            z = fPan.CG(i,2) - fPan.CG(j,2);
            coeffF = infcF(x, y, z, xC, yC, zC);
            mgAIC.CuYbwd(i,j) = coeffF[0];
            mgAIC.CvYbwd(i,j) = coeffF[1];
            mgAIC.CwYbwd(i,j) = coeffF[2];
        }
        for (int i = 0; i < fPan.nF; ++i) {
            x = fPan.CG(i,0) - fPan.CG(j,0);
            y = fPan.CG(i,1)+fPan.deltaMG - fPan.CG(j,1);
            z = fPan.CG(i,2) - fPan.CG(j,2);
            coeffF = infcF(x, y, z, xC, yC, zC);
            mgAIC.CuYfwd(i,j) = coeffF[0];
            mgAIC.CvYfwd(i,j) = coeffF[1];
            mgAIC.CwYfwd(i,j) = coeffF[2];
        }
        for (int i = 0; i < fPan.nF; ++i) {
            x = fPan.CG(i,0) - fPan.CG(j,0);
            y = fPan.CG(i,1) - fPan.CG(j,1);
            z = fPan.CG(i,2)-fPan.deltaMG - fPan.CG(j,2);
            coeffF = infcF(x, y, z, xC, yC, zC);
            mgAIC.CuZbwd(i,j) = coeffF[0];
            mgAIC.CvZbwd(i,j) = coeffF[1];
            mgAIC.CwZbwd(i,j) = coeffF[2];
        }
        for (int i = 0; i < fPan.nF; ++i) {
            x = fPan.CG(i,0) - fPan.CG(j,0);
            y = fPan.CG(i,1) - fPan.CG(j,1);
            z = fPan.CG(i,2)+fPan.deltaMG - fPan.CG(j,2);
            coeffF = infcF(x, y, z, xC, yC, zC);
            mgAIC.CuZfwd(i,j) = coeffF[0];
            mgAIC.CvZfwd(i,j) = coeffF[1];
            mgAIC.CwZfwd(i,j) = coeffF[2];
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
            x = fPan.CG(j,0) - fPan.CG(j,0);
            y = fPan.CG(j,1) + fPan.CG(j,1);
            z = fPan.CG(j,2) - fPan.CG(j,2);
            coeffF = infcF(x, y, z, xC, yC, zC);
            f2fAIC.Cu(j,j) += coeffF[0];
            f2fAIC.Cv(j,j) += coeffF[1];
            f2fAIC.Cw(j,j) += coeffF[2];
            for (int i = j+1; i < fPan.nF; ++i) {
                x = fPan.CG(i,0) - fPan.CG(j,0);
                y = fPan.CG(i,1) + fPan.CG(j,1);
                z = fPan.CG(i,2) - fPan.CG(j,2);
                coeffF = infcF(x, y, z, xC, yC, zC);
                f2fAIC.Cu(i,j) += coeffF[0];
                f2fAIC.Cv(i,j) += coeffF[1];
                f2fAIC.Cw(i,j) += coeffF[2];
                f2fAIC.Cu(j,i) += -coeffF[0];
                f2fAIC.Cv(j,i) += -coeffF[1];
                f2fAIC.Cw(j,i) += -coeffF[2];
            }

            // Field to field for minigrid
            for (int i = 0; i < fPan.nF; ++i) {
                x = fPan.CG(i,0)-fPan.deltaMG - fPan.CG(j,0);
                y = fPan.CG(i,1) + fPan.CG(j,1);
                z = fPan.CG(i,2) - fPan.CG(j,2);
                coeffF = infcF(x, y, z, xC, yC, zC);
                mgAIC.CuXbwd(i,j) += coeffF[0];
                mgAIC.CvXbwd(i,j) += coeffF[1];
                mgAIC.CwXbwd(i,j) += coeffF[2];
            }
            for (int i = 0; i < fPan.nF; ++i) {
                x = fPan.CG(i,0)+fPan.deltaMG - fPan.CG(j,0);
                y = fPan.CG(i,1) + fPan.CG(j,1);
                z = fPan.CG(i,2) - fPan.CG(j,2);
                coeffF = infcF(x, y, z, xC, yC, zC);
                mgAIC.CuXfwd(i,j) += coeffF[0];
                mgAIC.CvXfwd(i,j) += coeffF[1];
                mgAIC.CwXfwd(i,j) += coeffF[2];
            }
            for (int i = 0; i < fPan.nF; ++i) {
                x = fPan.CG(i,0) - fPan.CG(j,0);
                y = fPan.CG(i,1)-fPan.deltaMG + fPan.CG(j,1);
                z = fPan.CG(i,2) - fPan.CG(j,2);
                coeffF = infcF(x, y, z, xC, yC, zC);
                mgAIC.CuYbwd(i,j) += coeffF[0];
                mgAIC.CvYbwd(i,j) += coeffF[1];
                mgAIC.CwYbwd(i,j) += coeffF[2];
            }
            for (int i = 0; i < fPan.nF; ++i) {
                x = fPan.CG(i,0) - fPan.CG(j,0);
                y = fPan.CG(i,1)+fPan.deltaMG + fPan.CG(j,1);
                z = fPan.CG(i,2) - fPan.CG(j,2);
                coeffF = infcF(x, y, z, xC, yC, zC);
                mgAIC.CuYfwd(i,j) += coeffF[0];
                mgAIC.CvYfwd(i,j) += coeffF[1];
                mgAIC.CwYfwd(i,j) += coeffF[2];
            }
            for (int i = 0; i < fPan.nF; ++i) {
                x = fPan.CG(i,0) - fPan.CG(j,0);
                y = fPan.CG(i,1) + fPan.CG(j,1);
                z = fPan.CG(i,2)-fPan.deltaMG - fPan.CG(j,2);
                coeffF = infcF(x, y, z, xC, yC, zC);
                mgAIC.CuZbwd(i,j) += coeffF[0];
                mgAIC.CvZbwd(i,j) += coeffF[1];
                mgAIC.CwZbwd(i,j) += coeffF[2];
            }
            for (int i = 0; i < fPan.nF; ++i) {
                x = fPan.CG(i,0) - fPan.CG(j,0);
                y = fPan.CG(i,1) + fPan.CG(j,1);
                z = fPan.CG(i,2)+fPan.deltaMG - fPan.CG(j,2);
                coeffF = infcF(x, y, z, xC, yC, zC);
                mgAIC.CuZfwd(i,j) += coeffF[0];
                mgAIC.CvZfwd(i,j) += coeffF[1];
                mgAIC.CwZfwd(i,j) += coeffF[2];
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