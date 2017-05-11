//// Field structure
// Field is a container holding information about the field (i.e. the fluid domain enclosing the surface networks).
// It contains:
// - generic information (number of field panels,...)
// - geometrical information (field panels centers and vertices,...)
// - flow information (velocity,...)

#ifndef FPMV1_FIELD_H
#define FPMV1_FIELD_H

#include <Eigen/Dense>
#include <array>

struct Field {

	// Geometry
    int nX = 0; // Number of initial cells along X
	int nY = 0; // Number of initial cells along Y
	int nZ = 0; // Number of initial cells along Z
    int nF = 0; // Number of field cells
    Eigen::MatrixX3d CG; // Coordinates of field CGs
    Eigen::MatrixX2d vX, vY, vZ; // Coordinates of field vertices
    double deltaX, deltaY, deltaZ, deltaMG = 0; // Cell X, Y, Z and minigrid size
    int nE = 0; // Number of external field cells
	int nI = 0; // Number of internal field cells
	int nW = 0; // Number of field cells in the wake
    Eigen::VectorXi fMap; // Exterior/Interior map of the field
    Eigen::VectorXi eIdx, iIdx; // Exterior/Interior cell indices
    Eigen::VectorXi wMap; // Wake map of the field
    Eigen::VectorXi wIdx; // Wake cell indices

	// Flow
	Eigen::VectorXd sigma; // Field panel source singularity
	Eigen::MatrixX3d dSigma; // Field panel source singularity gradient
	Eigen::VectorXd sigmaTilda; // Modified field panel source singularity
	Eigen::MatrixX3d U; // Velocity components
	Eigen::VectorXd M; // Mach number
	Eigen::VectorXd rho; // Density
	Eigen::MatrixX3d dRho; // Density gradient
	Eigen::VectorXd a; // Speed of sound
};

#endif //FPMV1_FIELD_H
