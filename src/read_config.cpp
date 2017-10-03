//// Configuration reading
// Read data from /IO/*.cfg file (numerical parameters, symmetry, mach number, angle of attack, field cells) and put
// them into relevant variables
//
// I/O:
// - path: path to *.cfg file
// - numC: list of numerical parameters (structure)
// - symY: defines symmetry about Y axis
// - sRef: reference surface of the full wing
// - machInf: freestream Mach number
// - AoA: freestream angle of attack
// - box: temporary array to store vertices defining the domain
// - fPan: field panels (structure)

#include <iostream>
#include <fstream>
#include <sstream>
#include <array>
#include "read_config.h"

#define PI 3.14159

using namespace std;

void read_config(string path, Numerical_CST &numC, bool &symY, double &sRef, double &machInf, double &AoA, array<array<double, 3>, 8> &box,
                 Field &fPan, Subpanel &sp) {

    ifstream infile(path);
    string line;

    string name;
    bool symValue;
    double numValue;
    double geoValue;
    double flowValue;
    int gridValue;
    double ptsValue1, ptsValue2, ptsValue3;

    string symbol ; // " = "
    string comment;
    string comment_symbol; // " // "
    string section;

    bool secFLAG = 1;
    int i = 2;

    if(!infile.is_open())
    {
        cout << "File not found: " << path << endl;
        exit(EXIT_FAILURE);
    }

    getline(infile, line);
    cout << "Reading file... '" << line << "'"<< endl;

    while (getline(infile, line))
    {

        if (secFLAG) {
            stringstream ss(line);
            ss >> section;
            secFLAG = 0;
        }

        else {

            if (section == "$sym") {
                stringstream ss(line);
                ss >> name >> symbol >> symValue >> comment_symbol >> comment;

                if (name == "XZ_SYM") {
                    symY = symValue;
                    secFLAG = 1;
                }

                else {
                    cout << "Invalid parameter name: " << name << " at line " << i << endl;
                    exit(EXIT_FAILURE);
                }

            }

            else if (section == "$num") {
                stringstream ss(line);
                ss >> name >> symbol >> numValue >> comment_symbol >> comment;

                if (name == "boxTol") {
                    numC.TOLB = numValue;
                }
                else if (name == "resRed") {
                    numC.RRED = numValue;
                    secFLAG = 1;
                }
                else {
                    cout << "Invalid parameter name: " << name << " at line " << i << endl;
                    exit(EXIT_FAILURE);
                }
            }

            else if (section == "$geo") {
                stringstream ss(line);
                ss >> name >> symbol >> geoValue >> comment_symbol >> comment;

                if (name == "S_REF") {
                    sRef = geoValue;
                    secFLAG = 1;
                }

                else {
                    cout << "Invalid parameter name: " << name << " at line " << i << endl;
                    exit(EXIT_FAILURE);
                }
            }

            else if (section == "$flow") {
                stringstream ss(line);
                ss >> name >> symbol >> flowValue >> comment_symbol >> comment;

                if (name == "MACH") {
                    machInf = flowValue;
                }

                else if (name == "AoA") {
                    AoA = flowValue*PI/180;
                    secFLAG = 1;
                }
                else {
                    cout << "Invalid parameter name: " << name << " at line " << i << endl;
                    exit(EXIT_FAILURE);
                }
            }

            else if (section == "$grid") {
                stringstream ss(line);
                ss >> name >> symbol >> gridValue >> comment_symbol >> comment;

                if (name == "X_DIV") {
                    fPan.nX = gridValue;
                }

                else if (name == "Y_DIV") {
                    fPan.nY = gridValue;
                }

                else if (name == "Z_DIV") {
                    fPan.nZ = gridValue;
                    secFLAG = 1;
                }
                else {
                    cout << "Invalid parameter name: " << name << " at line " << i << endl;
                    exit(EXIT_FAILURE);
                }
            }

            else if (section == "$box") {
                stringstream ss(line);
                ss >> name >> symbol >> ptsValue1 >> ptsValue2 >> ptsValue3 >> comment_symbol >> comment;

                if (name == "BOX_X0Y0Z0") {
                    box[0][0] = ptsValue1; box[0][1] = ptsValue2; box[0][2] = ptsValue3;
                }

                else if (name == "BOX_X1Y0Z0") {
                    box[1][0] = ptsValue1; box[1][1] = ptsValue2; box[1][2] = ptsValue3;
                }

                else if (name == "BOX_X0Y0Z1") {
                    box[2][0] = ptsValue1; box[2][1] = ptsValue2; box[2][2] = ptsValue3;
                }

                else if (name == "BOX_X1Y0Z1") {
                    box[3][0] = ptsValue1; box[3][1] = ptsValue2; box[3][2] = ptsValue3;
                }

                else if (name == "BOX_X0Y1Z0") {
                    box[4][0] = ptsValue1; box[4][1] = ptsValue2; box[4][2] = ptsValue3;
                }

                else if (name == "BOX_X1Y1Z0") {
                    box[5][0] = ptsValue1; box[5][1] = ptsValue2; box[5][2] = ptsValue3;
                }

                else if (name == "BOX_X0Y1Z1") {
                    box[6][0] = ptsValue1; box[6][1] = ptsValue2; box[6][2] = ptsValue3;
                }

                else if (name == "BOX_X1Y1Z1") {
                    box[7][0] = ptsValue1; box[7][1] = ptsValue2; box[7][2] = ptsValue3;
                    secFLAG = 1;
                }
                else {
                    cout << "Invalid parameter name: " << name << " at line " << i << endl;
                    exit(EXIT_FAILURE);
                }
            }

            else if (section == "$sub") {
                stringstream ss(line);
                ss >> name >> symbol >> numValue >> comment_symbol >> comment;

                if (name == "NS") {
                    sp.NSs = (int) numValue;
                    sp.NS = sp.NSs * sp.NSs;
                }
                else if (name == "NC") {
                    sp.NC = numValue;
                }
                else if (name == "LC") {
                    sp.LC = numValue;
                    break;
                }
                else {
                    cout << "Invalid parameter name: " << name << " at line " << i << endl;
                    exit(EXIT_FAILURE);
                }
            }
            else {
                cout << "Invalid section name: " << section << " at line " << i-1 << endl;
                exit(EXIT_FAILURE);
            }
        }i++;
    }

    fPan.nF = fPan.nX * fPan.nY * fPan.nZ;

    cout << "Done reading config file!" << endl;
    cout << "XZ symmetry: " << symY << endl;
    cout << "Tolerance on box: " << numC.TOLB << endl;
    cout << "Order of magnitude of residual reduction: " << numC.RRED << endl;
    cout << "Reference surface: " << sRef << endl;
    cout << "Mach number: " << machInf << endl;
    cout << "Angle of attack: " << AoA*180/PI << endl;
    cout << "Number of cell (x): " << fPan.nX << endl;
    cout << "Number of cell (y): " << fPan.nY << endl;
    cout << "Number of cell (z): " << fPan.nZ << endl;
    cout << "Number of cells: " << fPan.nF << endl;
    cout << "Points defining the box:" << endl;
    for (int k = 0; k < 8; ++k)
        cout << box[k][0] << ' ' << box[k][1] << ' ' << box[k][2] << endl;
    cout << endl;
}