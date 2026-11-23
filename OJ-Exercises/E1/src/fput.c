#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "tools.h"
#include "fput.h"

/**
 * Calculate the acceleration (m = 1, kappa = 1, fixed ends u_0 = u_{N+1} = 0).
 *
 * Parameters
 * ----------
 *  accelerations - Vector where accelerations are to be written
 *  positions - Vector of positions
 *  alpha - Anharmonicity constant
 *  N - Number of atoms
 *
*/
void fput_acceleration(double *accelerations, double *positions, double alpha, const unsigned int N)
{
    for (unsigned int j = 0; j < N; ++j) {
        double u_left  = (j == 0)     ? 0.0 : positions[j - 1];
        double u_right = (j == N - 1) ? 0.0 : positions[j + 1];

        double d_right = u_right - positions[j];   /* extension of right spring */
        double d_left  = positions[j] - u_left;    /* extension of left spring  */

        accelerations[j] = (d_right - d_left) + alpha * (d_right * d_right - d_left * d_left);
    }
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
void fput_verlet_step(double *accelerations, double *positions, double *velocities,
                              double alpha, double timestep, const unsigned int N)
{
    /* v(t + dt/2) = v(t) + 0.5 a(t) dt ;  q(t + dt) = q(t) + v(t + dt/2) dt */
    for (unsigned int i = 0; i < N; ++i) {
        velocities[i] += 0.5 * accelerations[i] * timestep;
        positions[i]  += velocities[i] * timestep;
    }

    /* a(t + dt) from the new positions */
    fput_acceleration(accelerations, positions, alpha, N);

    /* v(t + dt) = v(t + dt/2) + 0.5 a(t + dt) dt */
    for (unsigned int i = 0; i < N; ++i) {
        velocities[i] += 0.5 * accelerations[i] * timestep;
    }
}

/* Energies, needed to check conservation (same units: m = kappa = 1) */
double fput_kinetic_energy(double *velocities, const unsigned int N)
{
    return 0.5 * dot_product(velocities, velocities, N);
}

double fput_potential_energy(double *positions, double alpha, const unsigned int N)
{
    double energy = 0.0;
    for (unsigned int i = 0; i <= N; ++i) {          /* N + 1 springs */
        double u_i  = (i == 0) ? 0.0 : positions[i - 1];
        double u_i1 = (i == N) ? 0.0 : positions[i];
        double d = u_i1 - u_i;
        energy += 0.5 * d * d + alpha * d * d * d / 3.0;
    }
    return energy;
}