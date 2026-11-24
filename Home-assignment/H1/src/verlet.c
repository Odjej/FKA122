
#include "potential.h"
/**
 * Calculate the acceleration.
 *
 *
 * Parameters
 * ----------
 *  accelerations - Vector where accelerations are to be written
 *  positions - Vector of positions
 *  alpha - Anharmonicity constant
 *  N - Number of atoms
 *
*/
void calculate_acceleration(double **accelerations, double **forces, const unsigned int N, double *masses)
{
    
    for(int i = 0; i < N; i++){
        for(int j = 0; j < 3; j++){
            accelerations[i][j] = forces[i][j] / masses[i];
        }
    }
}

double calculate_kinetic_energy(double **velocities, double *masses)
{
    const unsigned int N = 3;  // There are three particles
    double kinetic_energy = 0;
    for(int i = 0; i < N; i++){
        for(int j = 0; j < 3; j++){
        kinetic_energy += masses[i] * velocities[i][j] * velocities[i][j] / 2;
    }
    }
    // Return kinetic energy
    return kinetic_energy;
    


}

/**
 * Perform one velocity Verlet step
 *
 * Parameters
 * ----------
 *  accelerations - Vector of accelerations
 *  positions - Vector of positions
 *  velocities - Vector of velocities
 *  alpha - Anharmonicity constant
 *  timestep - Time step
 *  N - Number of atoms
 *
*/
void velocity_verlet_one_step(double **accelerations, double **forces, double **positions, double **velocities, double timestep, double lattice_param,const unsigned int N, double *masses)
{   
    for (int i = 0; i < N; i++) {
        for(int j = 0; j < 3; j++){
        velocities[i][j] += 0.5 * accelerations[i][j] * timestep;
        positions[i][j] += velocities[i][j] * timestep;
    }
    }
    
    get_forces_AL(forces, positions, N * lattice_param, N);
    calculate_acceleration(accelerations, forces, N, masses);
    
    for (int i = 0; i < N; i++) {
        for(int j = 0; j < 3; j++){
            velocities[i][j] += 0.5 * accelerations[i][j] * timestep;
    }
    }
    // Write to accelerations, positions and velocities vectors
}

