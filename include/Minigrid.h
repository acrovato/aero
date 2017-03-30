//// Minigrid structure
// Minigrid is a structure containing the flow variables at minigrid points

#ifndef FPMV1_MINIGRID_H
#define FPMV1_MINIGRID_H

#include <Eigen/Dense>

struct Minigrid {

    // Velocity
    Eigen::MatrixX3d UXbwd; // Velocity components @ X + deltaMG
    Eigen::MatrixX3d UXfwd; // Velocity components @ X - deltaMG
    Eigen::MatrixX3d UYbwd; // Velocity components @ Y + deltaMG
    Eigen::MatrixX3d UYfwd; // Velocity components @ Y - deltaMG
    Eigen::MatrixX3d UZbwd; // Velocity components @ Z + deltaMG
    Eigen::MatrixX3d UZfwd; // Velocity components @ Z - deltaMG

    // Density
    Eigen::VectorXd rhoXbwd; // Density @ X + deltaMG
    Eigen::VectorXd rhoXfwd; // Density @ X - deltaMG
    Eigen::VectorXd rhoYbwd; // Density @ Y + deltaMG
    Eigen::VectorXd rhoYfwd; // Density @ Y - deltaMG
    Eigen::VectorXd rhoZbwd; // Density @ Z + deltaMG
    Eigen::VectorXd rhoZfwd; // Density @ Z - deltaMG
};

#endif //FPMV1_MINIGRID_H
