#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "verlet.h"

/**
 * Calculate the acceleration.
 *
 * All vectors are assumed to be of length 3
 *
 * Parameters
 * ----------
 *  accelerations - Vector where accelerations are to be written
 *  positions - Vector of positions
 *  masses - Vector of masses
 *  kappa - Spring constant
 *
*/
void calculate_acceleration(double *accelerations, double *positions, double *masses, double kappa)
{
    // m1 q1'' = kappa * (q2 - q1)
    accelerations[0] = kappa * (positions[1] - positions[0]) / masses[0];
    // m2 q2'' = kappa * (q3 - 2 q2 + q1)
    accelerations[1] = kappa * (positions[2] - 2 * positions[1] + positions[0]) / masses[1];
    // m3 q3'' = kappa * (-q3 + q2)
    accelerations[2] = kappa * (-positions[2] + positions[1]) / masses[2];
}

/**
 * Calculate the potential energy
 *
 * All vectors are assumed to be of length 3
 *
 * Parameters
 * ----------
 *  positions - Vector of positions
 *  kappa - Spring constant
 * Returns
 * -------
 *  Potential energy
 *
*/
double calculate_potential_energy(double *positions, double kappa)
{
    const unsigned int N = 3;  // There are three particles
    double potential_energy = 0;
    for(int i = 0; i < N - 1; i++){
        potential_energy += 0.5 * kappa * (positions[i+1] - positions[i]) * (positions[i+1] - positions[i]);
    }
    return potential_energy;
}

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
double calculate_kinetic_energy(double *velocities, double *masses)
{
    const unsigned int N = 3;  // There are three particles

    double kinetic_energy = 0;
    for(int i = 0; i < N; i++){
        kinetic_energy += 0.5 * masses[i] * velocities[i] * velocities[i];
    }
    return kinetic_energy;
}

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
void velocity_verlet_one_step(double *accelerations, double *positions, double *velocities,
                              double *masses, double kappa, double timestep)
{
    const unsigned int N = 3;  // There are three particles
    // v_i(t + dt/2) = v_i(t) + 0.5 * a_i(t) * dt
    // q_i(t + dt)   = q_i(t) + v_i(t + dt/2) * dt
    for (int i = 0; i < N; i++) {
        velocities[i] += 0.5 * accelerations[i] * timestep;
        positions[i] += velocities[i] * timestep;
    }

    // Calculate new accelerations a_i(t + dt) from the new positions
    calculate_acceleration(accelerations, positions, masses, kappa);

    // v_i(t + dt) = v_i(t + dt/2) + 0.5 * a_i(t + dt) * dt
    for (unsigned int i = 0; i < N; i++) {
        velocities[i] += 0.5 * accelerations[i] * timestep;
    }
}
