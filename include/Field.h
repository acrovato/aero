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
    int nX; // Number of initial cells along X
	int nY; // Number of initial cells along Y
	int nZ; // Number of initial cells along Z
    int nF; // Number of field cells
    Eigen::MatrixX3d CG; // Coordinates of field CGs
    Eigen::MatrixX2d vX, vY, vZ; // Coordinates of field vertices
    double deltaX, deltaY, deltaZ, deltaMG; // Cell X, Y, Z and minigrid size
    int nE; // Number of external field cells
	int nI; // Number of internal field cells
    Eigen::VectorXi fMap; // Exterior/Interior map of the field
    Eigen::VectorXi eIdx, iIdx; // Exterior/Interior cell indices
    Eigen::VectorXi wMap; // Exterior/Interior map of the field
	
	// Flow
	Eigen::VectorXd sigma; // Field panel source singularity
	Eigen::MatrixX3d U; // Velocity components
	Eigen::VectorXd M; // Mach number
	Eigen::VectorXd rho; // Density
	Eigen::VectorXd a; // Speed of sound
};

#endif //FPMV1_FIELD_H
