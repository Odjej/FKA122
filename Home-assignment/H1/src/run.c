#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
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

double calculate_pressure(double volume, double kinetic_energy, double virial)
{
    return (2.0 * kinetic_energy / 3.0 + virial) / volume;
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



static int run_default(void)
{
    int N = 4;
    double lattice_param = 4.03230237;
    int n_atoms = 4 * N * N * N;
    int total_timesteps = 200000;
    double dt = .002;
    
    double T_eq = 1200 + 273.15;
    double T_aim = 700 + 273.15;
    bool melt = true;


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
    double virial;
    
    bool steady_state = false;

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


    calculate(&potential_energy, &virial, forces, positions, N * lattice_param, n_atoms);
    calculate_acceleration(accelerations, forces, n_atoms, masses);
    
    FILE *file = fopen("data/pressure_temperature_over_time_fixed units.csv", "w");
    
    FILE *files[10];
    char filename[50];

    for(int i = 0; i < 10; i++){
        snprintf(filename, sizeof(filename), "data/positions_%d.csv", i);
        files[i] = fopen(filename, "w");
        fprintf(files[i], "time,x,y,z,steady_state\n");
    }
    
    fprintf(file, "time,potential_energy,kinetic_energy,total_energy,temperature,pressure,lattice_param,steady_state\n");

    printf("Running simulations \n");
    for (int i = 0; i < total_timesteps; i++) {
        

        velocity_verlet_one_step(accelerations, forces, positions, velocities, dt, lattice_param * N, n_atoms, masses, &potential_energy, &virial);
        kinetic_energy = calculate_kinetic_energy(velocities, masses, n_atoms);
        total_energy = potential_energy + kinetic_energy;
        

        pressure = calculate_pressure(pow(N * lattice_param, 3), kinetic_energy, virial);
        temperature = calculate_temperature(n_atoms, velocities, masses);
        
        if(steady_state == false){
            if(i > 100000){
                steady_state = true;
            }}

        if(steady_state == false){
        scale_positions(positions, n_atoms, dt, tau_P, pressure, P_eq, kappa_T);
        scale_speeds(velocities, n_atoms, dt, tau_T, temperature, T_eq);
        scale_lattice_param(&lattice_param, dt, tau_P, pressure, P_eq, kappa_T);
        }

        calculate(&potential_energy, NULL, forces, positions, N * lattice_param, n_atoms);
        calculate_acceleration(accelerations, forces, n_atoms, masses);

        fprintf(file, "%.4f,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f,%d\n", dt * (1 + i) ,potential_energy, kinetic_energy, total_energy, temperature, pressure * pressure_conversion / 1e9, lattice_param, steady_state);
        if ((i + 1) % 2000 == 0){
            if(T_eq > T_aim && melt){
                T_eq = fmax(T_aim, T_eq - 50);
            }
            printf("Timestep %d done \n", (1 + i));

            for(int j = 0; j < 10; j++){
                fprintf(files[j], "%.4f,%.12f,%.12f,%.12f,%d\n", dt * (1 + i), positions[j*10][0], positions[j*10][1], positions[j*10][2], steady_state);
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

static int run_task1(void)
{
    int N = 4;
    int n_atoms = 4 * N * N * N;
    double **positions = create_2D_array(n_atoms, 3);
    FILE *file = fopen("data/lattice_potential_energy.csv", "w");

    for (int i = 1; i < 200; i++) {
        double lattice_param = 3.5 + 0.01 * i;
        init_fcc(positions, N, lattice_param);

        double potential_energy = get_energy_AL(positions, N * lattice_param, n_atoms);
        fprintf(file, "%.6f,%.6f,\n", lattice_param, potential_energy);
    }

    fclose(file);
    destroy_2D_array(positions, n_atoms);
    return 0;
}

static int run_task2(void)
{
    int N = 4;
    double lattice_param = 4.03230237;
    int n_atoms = 4 * N * N * N;
    const double AL_mass = 0.00279640;
    const double fs = 1e-3;   // 1 fs = 0.001 ps

    double dt_fs[] = {0.5, 1, 2, 4, 6, 8, 12, 16, 18, 20, 24};
    int n_dt = sizeof(dt_fs) / sizeof(dt_fs[0]);
    double t_total = 5000.0 * fs;        // same physical time (5 ps) for every dt

    double **positions = create_2D_array(n_atoms, 3);
    double **forces = create_2D_array(n_atoms, 3);
    double **accelerations = create_2D_array(n_atoms, 3);
    double **velocities = create_2D_array(n_atoms, 3);
    double *masses = malloc(sizeof(double) * n_atoms);


    gsl_rng_env_setup();
    gsl_rng *r = gsl_rng_alloc(gsl_rng_default);

    for (int k = 0; k < n_dt; k++) {
        double dt = dt_fs[k] * fs;
        int total_timesteps = (int)(t_total / dt);

        // Identical starting configuration for every dt
        init_fcc(positions, N, lattice_param);
        gsl_rng_set(r, 42);
        for (int i = 0; i < n_atoms; i++) {
            masses[i] = AL_mass;
            for (int j = 0; j < 3; j++) {
                positions[i][j] += lattice_param * 0.065 * (1 - 2 * gsl_rng_uniform(r));
                velocities[i][j] = 0;
            }
        }

        double potential_energy, kinetic_energy, total_energy;
        calculate(&potential_energy, NULL, forces, positions,
                  N * lattice_param, n_atoms);
        calculate_acceleration(accelerations, forces, n_atoms, masses);

        char filename[128];
        snprintf(filename, sizeof(filename), "data/energies_dt_%.1ffs.csv", dt_fs[k]);
        FILE *file = fopen(filename, "w");
        if (!file) { perror("fopen"); return 1; }
        fprintf(file, "time_fs,potential_energy,kinetic_energy,total_energy\n");

        // Log t = 0 so we can plot E_tot(t) - E_tot(0)
        kinetic_energy = calculate_kinetic_energy(velocities, masses, n_atoms);
        fprintf(file, "%.6e,%.12f,%.12f,%.12f\n", 0.0, potential_energy,
                kinetic_energy, potential_energy + kinetic_energy);

        printf("dt = %.1f fs, %d steps\n", dt_fs[k], total_timesteps);
        for (int i = 0; i < total_timesteps; i++) {
            velocity_verlet_one_step(accelerations, forces, positions, velocities, dt,
                                     lattice_param * N, n_atoms, masses,
                                     &potential_energy, NULL);
            kinetic_energy = calculate_kinetic_energy(velocities, masses, n_atoms);
            total_energy = potential_energy + kinetic_energy;

            fprintf(file, "%.6e,%.12f,%.12f,%.12f\n", dt_fs[k] * (i + 1),
                    potential_energy, kinetic_energy, total_energy);

            // Stop early if the run has blown up
            if (!isfinite(total_energy)) {
                printf("  energy diverged at step %d, stopping this dt\n", i + 1);
                break;
            }
        }
        fclose(file);
    }

    destroy_2D_array(positions, n_atoms);
    destroy_2D_array(velocities, n_atoms);
    destroy_2D_array(forces, n_atoms);
    destroy_2D_array(accelerations, n_atoms);
    free(masses);
    gsl_rng_free(r);
    return 0;
}

int run(int argc, char *argv[])
{
    if (argc == 1) {
        return run_default();
    }
    if (argc == 2) {
        if (strcmp(argv[1], "task1") == 0) {
            return run_task1();
        }
        if (strcmp(argv[1], "task2") == 0) {
            return run_task2();
        }
    }

    fprintf(stderr, "Usage: %s [task1|task2]\n", argv[0]);
    return EXIT_FAILURE;
}
