#pragma once

/**
 * Calculate the acceleration.
 *
 * All matrices are assumed to be of n_atoms x 3
 *
 * Parameters
 * ----------
 *  accelerations - Vector where accelerations are to be written
 *  forces - Vector of forces
 *  N - number of atoms
 *  masses - Vector of masses
 *
*/
void 
calculate_acceleration(double **accelerations, double **forces, const unsigned int N, double *masses);



/**
 * Calculate the kinetic energy
 *
 * All vectors are assumed to be of length 3
 *
 * Parameters
 * ----------
 *  velocities - Vector of velocities
 *  masses - Vector of masses
 *
 * Returns
 * -------
 *  Kinetic energy
 *
*/
double 
calculate_kinetic_energy(double **velocities, double *masses);


/**
 * Perform one velocity Verlet step
 *
 * All vectors are assumed to be of length 3
 *
 * Parameters
 * ----------
 *  accelerations - Vector of accelerations
 *  positions - Vector of positions
 *  velocities - Vector of velocities
 *  masses - Vector of masses
 *  kappa - Spring constant
 *  timestep - Time step
 *
*/
void 
velocity_verlet_one_step(double **accelerations, double **forces, double **positions, double **velocities, double timestep, double lattice_param,const unsigned int N, double *masses);
