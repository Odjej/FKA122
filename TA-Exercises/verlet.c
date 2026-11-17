void calculate_acceleration(double *accelerations, double *positions, double *masses, double kappa)
{
    const unsigned int N = 3;  // There are three particles
    // Write to the accelerations vector
    accelerations[0] = kappa * (positions[1] - positions[0]) / masses[0];
    accelerations[1] = kappa * (positions[2] - 2 * positions[1] + positions[0])/ masses[1];
    accelerations[2] = kappa * (-positions[2] + positions[1]) / masses[2];
}

double calculate_potential_energy(double *positions, double kappa)
{
    const unsigned int N = 3;  // There are three particles
    double potential_energy = 0;
    for(int i = 0; i < N - 1; i++){
        potential_energy += (kappa / 2) * (positions[i+1] - positions[i]) * (positions[i+1] - positions[i]);
    }

    // Return potential energy
    return potential_energy;
}

double calculate_kinetic_energy(double *velocities, double *masses)
{
    const unsigned int N = 3;  // There are three particles
    double kinetic_energy = 0;
    for(int i = 0; i < N; i++){
        kinetic_energy += masses[i] * velocities[i] * velocities[i] / 2;
    // Return kinetic energy
}
    return kinetic_energy;

}

void velocity_verlet_one_step(double *accelerations, double *positions, double *velocities,
                              double *masses, double kappa, double timestep)
{
    const unsigned int N = 3;  // There are three particles
    for (int i = 0; i < N; i++) {
        velocities[i] += 0.5 * accelerations[i] * timestep;
        positions[i] += velocities[i] * timestep;
    }
    
    calculate_acceleration(accelerations, positions, masses, kappa);

    for (int i = 0; i < N; i++) {
        velocities[i] += 0.5 * accelerations[i] * timestep;
    }
    // Write to accelerations, positions and velocities vectors
}