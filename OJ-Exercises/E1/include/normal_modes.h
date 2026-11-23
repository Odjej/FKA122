#pragma once
void transform_to_normal_modes(double *positions, double *Q, const unsigned int N);
void calculate_normal_mode_energies(double *energies, double *positions, double *velocities, const unsigned int N);