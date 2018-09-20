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

//// Network structure
// Network is a container holding information about a network (i.e. a group of surface panels).
// It contains:
// - generic information (number of panels,...)
// - geometrical information (Panels centers and vertices,...)
// - flow information ( velocity,...)

#ifndef FPMV1_NETWORK_H
#define FPMV1_NETWORK_H

#include <Eigen/Dense>

struct Network {

	// Geometry 
    int nC = 0; // Number of chordwise points
    int nS = 0; // Number of spanwise points
	int nC_ = 0; // Number of chordwise panel
	int nS_ = 0; // Number of spanwise panel
    int nP = 0; // Number of panels
    Eigen::MatrixX3d CG; // Coordinates of panel CGs
    Eigen::MatrixX3d v0; // Coordinates of panel first vertex
    Eigen::MatrixX3d v1; // Coordinates of panel second vertex
    Eigen::MatrixX3d v2; // Coordinates of panel third vertex
    Eigen::MatrixX3d v3; // Coordinates of panel fourth vertex
    Eigen::MatrixX3d l; // Components of panel (first) longitudinal vector
    Eigen::MatrixX3d t; // Components of panel (second) longitudinal vector
    Eigen::MatrixX3d p; // Components of panel perpendicular vector
    Eigen::MatrixX3d n; // Components of panel normal vector
    Eigen::VectorXd S; // Panel surface

	// Flow
	Eigen::VectorXd tau; // Surface panel source singularity
	Eigen::VectorXd mu; // Surface panel doublet singularity
	Eigen::MatrixX3d U; // Velocity components
	Eigen::VectorXd M; // Mach number
	Eigen::VectorXd cP; //Pressure coefficient
	};

#endif //FPMV1_NETWORK_H
