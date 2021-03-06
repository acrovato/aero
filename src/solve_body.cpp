/* 
// Copyright 2018 University of Liege
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// 
//     http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
// Authors:
// - Adrien Crovato
*/

//// Incompressible module (body solver)
// Compute sources on the body surface from the updated  impermeability B.C. and solve for doublets on the body surface
//
// Ref: Katz & Plotkin (2001), Low-Speed Aerodynamics
// 1) tau = - (U_inf + Usigma) \cdot n
// 2) mu = A \ (B * tau)
//
// I/O:
// - vInf: freestream velocity vector
// - RHS: right-hand side of the linear system of equations
// - vSigma: field source induced body velocity
// - bPan: (network of) body panels (structure)
// - b2bAIC: body to body AIC (structure)

#include <iostream>
#include <Eigen/Dense>
#include "solve_body.h"

using namespace std;
using namespace Eigen;

void solve_body(Vector3d &vInf, VectorXd &RHS, MatrixX3d &vSigma, Network &bPan, Body_AIC &b2bAIC) {

    //// Begin
    cout << "Computing surface singularities... " << flush;

    //// Update B.C.
    // BC: tau_i = n_i * V_inf + n_i * V_sigma
    for (int i = 0; i < bPan.nP; ++i)
        bPan.tau(i) = - bPan.n.row(i).dot(vInf.transpose()) - bPan.n.row(i).dot(vSigma.row(i));

    //// Solve
    // Compute RHS
    RHS = - b2bAIC.B * bPan.tau;
    // Solve: A*mu + B*sigma = 0 (Dirichlet: interior potential = 0)
    #ifdef ON_UNIX
        bPan.mu = b2bAIC.A.householderQr().solve(RHS);
    #else
        bPan.mu = b2bAIC.A.fullPivLu().solve(RHS);
    #endif

    //// Control display
    cout << "Done!" << endl;
    cout << "Surface sources strength: " << bPan.tau.norm() << endl;
    #ifdef VERBOSE
        cout << "Sources min. value: " << bPan.tau.minCoeff() << endl;
        cout << "Sources max. value: " << bPan.tau.maxCoeff() << endl;
        for (int i = 0; i < bPan.nP; ++i)
            cout << i << ' ' << bPan.tau(i) << endl;
    #endif
    cout << "Surface doublets strength: " << bPan.mu.norm() << endl;
    #ifdef VERBOSE
        cout << "Doublets min. value: " << bPan.mu.minCoeff() << endl;
        cout << "Doublets max. value: " << bPan.mu.maxCoeff() << endl;
        for (int i = 0; i < bPan.nP; ++i)
            cout << i << ' ' << bPan.mu(i) << endl;
            cout << "System solved with relative error: " << (b2bAIC.A * bPan.mu - RHS).norm() / RHS.norm() << endl;
    #endif
    cout << endl;
}
