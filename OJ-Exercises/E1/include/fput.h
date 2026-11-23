#pragma once
#include <stdio.h>
void fput_acceleration(double *accelerations, double *positions, double alpha, const unsigned int N);

void fput_verlet_step(double *accelerations, double *positions, double *velocities,
                      double alpha, double timestep, const unsigned int N);

double fput_kinetic_energy(double *velocities, const unsigned int N);

double fput_potential_energy(double *positions, double alpha, const unsigned int N);