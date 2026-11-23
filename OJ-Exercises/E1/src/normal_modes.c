#include <stdlib.h>
#include <math.h>

#include "tools.h"
#include "normal_modes.h"

/* Build the N x N transform matrix S_ki = sqrt(2/(N+1)) sin(i k pi/(N+1)), 1-indexed in the formula */
static double **create_transform_matrix(const unsigned int N)
{
    double **S = create_2D_array(N, N);
    if (S == NULL)
        return NULL;

    const double prefactor = sqrt(2.0 / (N + 1));
    for (unsigned int k = 0; k < N; ++k)
        for (unsigned int i = 0; i < N; ++i)
            S[k][i] = prefactor * sin((i + 1) * (k + 1) * M_PI / (N + 1));
    return S;
}

/**
 * Transform to normal modes (m = 1). Also works as the inverse transform,
 * since the matrix is symmetric and orthogonal.
 */
void transform_to_normal_modes(double *positions, double *Q, const unsigned int N)
{
    double **S = create_transform_matrix(N);
    if (S == NULL)
        return;

    matrix_vector_multiplication(Q, S, positions, N, N);
    destroy_2D_array(S, N);
}

/**
 * Calculate the normal mode energies E_k = 0.5 (P_k^2 + omega_k^2 Q_k^2)
 */
void calculate_normal_mode_energies(double *energies, double *positions, double *velocities, const unsigned int N)
{
    double *Q = malloc(sizeof(double) * N);
    double *P = malloc(sizeof(double) * N);
    if (Q == NULL || P == NULL) {
        free(Q);
        free(P);
        return;
    }

    transform_to_normal_modes(positions, Q, N);
    transform_to_normal_modes(velocities, P, N);

    for (unsigned int k = 0; k < N; ++k) {
        double omega = 2.0 * sin((k + 1) * M_PI / (2.0 * (N + 1)));   /* kappa = m = 1 */
        energies[k] = 0.5 * (P[k] * P[k] + omega * omega * Q[k] * Q[k]);
    }

    free(Q);
    free(P);
}