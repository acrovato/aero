//// Field to field influence coefficient computation
// Compute field source potential influence coefficient between 2 field panels
//
// Reference:
// Zakir F. Seidov, P.I. Skvirsky (2000), Gravitational potential and energy of homogeneous rectangular parallelepiped,
// Department of Physics, Ben-Gurion University of the Negev, Beer-Sheva, 84105, Israel
//
// Inputs:
// - x, y, z: coordinates of target point
// - xC: x coordinates of influencing field cell vertices
// - yC: y coordinates of influencing field cell vertices
// - zC: z coordinates of influencing field cell vertices
//
// Output:
// - coeff: AIC

#include <iostream>
#include <array>
#include <cmath>
#include "infcFF.h"
#define PI 3.14159
using namespace std;

double infcFF(double x, double y, double z, array<double,2> &xC, array<double,2> &yC, array<double,2> &zC) {

    double coeff = 0;
    double A, B, C, I, II, III, R;

    for (int i = 1; i <= 2; ++i) {
        for (int j = 1; j <= 2; ++j) {
            for (int k = 1; k <= 2; ++k) {
                A = x - xC[i-1];
                B = y - yC[j-1];
                C = z - zC[k-1];
                R = sqrt(A*A + B*B + C*C);
                I = B*C*log(A+R) - A*A/2 * atan((B*C)/(A*R));
                II = C*A*log(B+R) - B*B/2 * atan((C*A)/(B*R));
                III = A*B*log(C+R) - C*C/2 * atan((A*B)/(C*R));
                coeff -= pow((-1), (double)(i+j+k)) * (I+II+III);
            }
        }
    }
    coeff *= -1/(4*PI);

    return coeff;
}