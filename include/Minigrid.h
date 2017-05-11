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

    Eigen::MatrixX3d UXmidbwd; // Velocity components @ X + deltaMG/2
    Eigen::MatrixX3d UXmidfwd; // Velocity components @ X - deltaMG/2
    Eigen::MatrixX3d UYmidbwd; // Velocity components @ Y + deltaMG/2
    Eigen::MatrixX3d UYmidfwd; // Velocity components @ Y - deltaMG/2
    Eigen::MatrixX3d UZmidbwd; // Velocity components @ Z + deltaMG/2
    Eigen::MatrixX3d UZmidfwd; // Velocity components @ Z - deltaMG/2

    // Density
    Eigen::VectorXd rhoXbwd; // Density @ X + deltaMG
    Eigen::VectorXd rhoXfwd; // Density @ X - deltaMG
    Eigen::VectorXd rhoYbwd; // Density @ Y + deltaMG
    Eigen::VectorXd rhoYfwd; // Density @ Y - deltaMG
    Eigen::VectorXd rhoZbwd; // Density @ Z + deltaMG
    Eigen::VectorXd rhoZfwd; // Density @ Z - deltaMG

    Eigen::VectorXd rhoXmidbwd; // Density @ X + deltaMG/2
    Eigen::VectorXd rhoXmidfwd; // Density @ X - deltaMG/2
    Eigen::VectorXd rhoYmidbwd; // Density @ Y + deltaMG/2
    Eigen::VectorXd rhoYmidfwd; // Density @ Y - deltaMG/2
    Eigen::VectorXd rhoZmidbwd; // Density @ Z + deltaMG/2
    Eigen::VectorXd rhoZmidfwd; // Density @ Z - deltaMG/2

    // Density gradient
    Eigen::MatrixX3d dRhoXmidbwd; // Density gradient @ X + deltaMG/2
    Eigen::MatrixX3d dRhoXmidfwd; // Density gradient @ X - deltaMG/2
    Eigen::MatrixX3d dRhoYmidbwd; // Density gradient @ Y + deltaMG/2
    Eigen::MatrixX3d dRhoYmidfwd; // Density gradient @ Y - deltaMG/2
    Eigen::MatrixX3d dRhoZmidbwd; // Density gradient @ Z + deltaMG/2
    Eigen::MatrixX3d dRhoZmidfwd; // Density gradient @ Z - deltaMG/2

    // Field source
    Eigen::VectorXd sigmaXmidbwd; // Field source @ X + deltaMG/2
    Eigen::VectorXd sigmaXmidfwd; // Field source @ X - deltaMG/2
    Eigen::VectorXd sigmaYmidbwd; // Field source @ Y + deltaMG/2
    Eigen::VectorXd sigmaYmidfwd; // Field source @ Y - deltaMG/2
    Eigen::VectorXd sigmaZmidbwd; // Field source @ Z + deltaMG/2
    Eigen::VectorXd sigmaZmidfwd; // Field source @ Z - deltaMG/2
};

#endif //FPMV1_MINIGRID_H
