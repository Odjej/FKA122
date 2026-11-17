#include "verlet.h"
#include "tools.h"
#include "unit_conversion.h"
#include <stdlib.h>
#include <stdio.h>

void save_positions(char *filename, double **matrix, int rows)
    {
        FILE *file = fopen(filename, "w");

        for (int i = 0; i < rows; i++){
            for (int j = 0; j < 3; j++){
                fprintf(file, "%.12f", matrix[i][j]);

                if(j < 2){
                    fprintf(file, ",");
                }
            }
            fprintf(file, "\n");
        }
        fclose(file);
    }

    
void save_energies(char *filename, double *matrix, int rows)
    {
        FILE *file = fopen(filename, "w");

        for (int i = 0; i < rows; i++){
                fprintf(file, "%.12e", matrix[i]);

            
            fprintf(file, "\n");
        }
        fclose(file);
    }

int
main(void){
    int total_timesteps = 1e5;
    int N = 3;
    double timestep = 0.0001;

    double spring_constant = 1600.0;
    double kappa = spring_constant_to_simulation_units(spring_constant);

    double *accelerations = malloc(N * sizeof(double));
    double *velocities = malloc(N * sizeof(double));
    double *masses = malloc(N * sizeof(double));

    double **positions = create_2D_array(total_timesteps, N);
    double *potential_energy = malloc(total_timesteps * sizeof(double));
    double *kinetic_energy = malloc(total_timesteps * sizeof(double));

    positions[0][0] = 0.01;
    positions[0][1] = 0.005;
    positions[0][2] = -0.005;

    velocities[0] = 0;
    velocities[1] = 0;
    velocities[2] = 0;
    
    double carbon_mass_amu = 12.011;
    double oxygen_mass_amu = 15.999;

    double carbon_mass = amu_to_simulation_units(carbon_mass_amu);
    double oxygen_mass = amu_to_simulation_units(oxygen_mass_amu);

    masses[0] = oxygen_mass;
    masses[1] = carbon_mass;
    masses[2] = oxygen_mass;
    
    calculate_acceleration(accelerations, positions[0], masses, kappa);

    for (int i = 0; i < total_timesteps; i++) {

        potential_energy[i] = calculate_potential_energy(positions[i], kappa);
        kinetic_energy[i] = calculate_kinetic_energy(velocities, masses);

        if (i < total_timesteps - 1) {

            for (int j = 0; j < N; j++) {
                positions[i + 1][j] = positions[i][j];
            }

            velocity_verlet_one_step(
                accelerations, positions[i + 1],
                velocities, masses, kappa, timestep
            );
        }
    }
    save_energies("potential_energy.csv", potential_energy, total_timesteps);
    save_energies("kinetic_energy.csv", kinetic_energy, total_timesteps);
    
    save_positions("positions.csv", positions, total_timesteps);
    
    
    double *spectrum = malloc(total_timesteps * sizeof(double));
    double *frequencies = malloc(total_timesteps * sizeof(double));
    
    powerspectrum(spectrum, potential_energy, total_timesteps, timestep);
    fft_freq(frequencies, total_timesteps, timestep);

    save_energies("powerspectrum.csv", spectrum, total_timesteps);
    save_energies("fft_frqs.csv", frequencies, total_timesteps);
    destroy_2D_array(positions, total_timesteps);

    free(accelerations);
    free(velocities);
    free(masses);
    free(potential_energy);
    free(kinetic_energy);
    free(spectrum);
    free(frequencies);

    
    return 0;
}