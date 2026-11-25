#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "lattice.h"
#include "potential.h"
#include "tools.h"
#include "verlet.h"

void 
scale_speeds(double **velocities, int n_atoms, double dt, double tau_T, double T_t, double T_eq ){
    double alpha_T = 1. + 2.*dt / (tau_T) * (T_eq - T_t) / T_t;

    for(int i = 0; i < n_atoms; i++){
        for(int j = 0; j < 3; j++){
            velocities[i][j] *= sqrt(alpha_T);
        }
    }
}

void 
scale_positions(double **positions, int n_atoms, double dt, double tau_P, double P_t, double P_eq, double kappa_T){
    double alpha_P = 1 - kappa_T*dt / (tau_P) * (P_eq - P_t);

    for(int i = 0; i < n_atoms; i++){
        for(int j = 0; j < 3; j++){
            positions[i][j] *= cbrt(alpha_P); // cubic root
        }
    }
}

void
scale_lattice_param(double *lattice_param, double dt, double tau_P, double P_t, double P_eq, double kappa_T){
    double alpha_P = 1 - kappa_T*dt / (tau_P) * (P_eq - P_t);
    *lattice_param *= cbrt(alpha_P);
}

double calculate_pressure(double volume, double n_atoms, double **forces, double **positions, double **velocities, double *masses){
    double P_sum = 0.0;
    for(int i = 0; i < n_atoms; i++){
        for(int j = 0; j < 3; j++){
        P_sum += masses[i] * velocities[i][j] *  velocities[i][j] + positions[i][j] * forces[i][j];
    }
}
    return P_sum / (3.0 * volume);
}

double calculate_temperature(double n_atoms, double **velocities, double *masses){
    double T_sum = 0.0;
    double boltzmann = 8.617333262E-5;
    for(int i = 0; i < n_atoms; i++){
        for(int j = 0; j < 3; j++){
            T_sum += masses[i] * velocities[i][j] *  velocities[i][j];
        }
    }
    return T_sum / (3 * n_atoms * boltzmann);
}



int
run(
    int argc,
    char *argv[]
   )
{
    int N = 4;
    double lattice_param = 4.03;
    int n_atoms = 4 * N * N * N;
    int total_timesteps = 10000;
    double dt = .005;
    
    double T_eq = 500 + 273.15;
    
    double temperature;
    double pressure;

    double tau_T = dt * 300;
    double tau_P = dt * 300;

    double pressure_conversion = 1.602176634E11;

    double P_eq = 1.E5 / pressure_conversion;
    double kappa_T = pressure_conversion / 76.E9;

    double **positions = create_2D_array(n_atoms, 3);
    double **forces = create_2D_array(n_atoms, 3);
    double **accelerations = create_2D_array(n_atoms, 3);
    double **velocities = create_2D_array(n_atoms, 3);
    
    double *masses = malloc(sizeof(double) * n_atoms);
    
    double potential_energy;
    double kinetic_energy;
    double total_energy;

    const double AL_mass = 0.00279640;

    
    init_fcc(positions, N, lattice_param);
    
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    int seed = 42;
    gsl_rng_set(r, seed);

    for(int i = 0; i < n_atoms; i++){
        masses[i] = AL_mass;
        for(int j = 0; j < 3; j++){
            positions[i][j] += lattice_param * 0.065 * (1- 2 *gsl_rng_uniform(r));
            velocities[i][j] = 0;
        }
    }


    get_forces_AL(forces, positions, N * lattice_param, n_atoms);
    calculate_acceleration(accelerations, forces, n_atoms, masses);
    
    FILE *file = fopen("pressure_temperature_over_time_fixed units.csv", "w");
    
    FILE *files[10];
    char filename[50];

    for(int i = 0; i < 10; i++){
        snprintf(filename, sizeof(filename), "positions_%d.csv", i);
        files[i] = fopen(filename, "w");
        fprintf(files[i], "time,x,y,z\n");
    }
    
    fprintf(file, "time,potential_energy,kinetic_energy,total_energy,temperature,pressure,lattice_param\n");

    printf("Running simulations \n");
    for (int i = 0; i < total_timesteps; i++) {
        
        velocity_verlet_one_step(accelerations, forces, positions, velocities, dt, lattice_param * N, n_atoms, masses);
        potential_energy = get_energy_AL(positions, N * lattice_param, n_atoms);
        kinetic_energy = calculate_kinetic_energy(velocities, masses, n_atoms);
        total_energy = potential_energy + kinetic_energy;
        
        pressure = calculate_pressure((N * lattice_param) * (N * lattice_param) * (N * lattice_param), n_atoms, forces, positions, velocities, masses);
        temperature = calculate_temperature(n_atoms, velocities, masses);

        scale_positions(positions, n_atoms, dt, tau_P, pressure, P_eq, kappa_T);
        scale_speeds(velocities, n_atoms, dt, tau_T, temperature, T_eq);
        scale_lattice_param(&lattice_param, dt, tau_P, pressure, P_eq, kappa_T);
        fprintf(file, "%.4f,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f \n", dt * (1 + i) ,potential_energy, kinetic_energy, total_energy, temperature, pressure * pressure_conversion / 1e9, lattice_param);
        
        if ((i + 1) % 100 == 0){
            printf("Timestep %d done \n", (1 + i));
            for(int j = 0; j < 10; j++){
                fprintf(files[j], "%.4f,%.12f,%.12f,%.12f\n", dt * (1 + i), positions[j*10][0], positions[j*10][1], positions[j*10][2]);
            }
        }
    }
    for(int i = 0; i < 10; i++){
        fclose(files[i]);
        }
    fclose(file);

    destroy_2D_array(positions, n_atoms);
    destroy_2D_array(velocities, n_atoms);
    destroy_2D_array(forces, n_atoms);
    destroy_2D_array(accelerations, n_atoms);
    free(masses);
    gsl_rng_free(r);
    return 0;
    }
