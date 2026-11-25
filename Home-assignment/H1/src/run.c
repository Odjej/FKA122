#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "lattice.h"
#include "potential.h"
#include "tools.h"
#include "verlet.h"

double 
scale_speeds(double **velocities, int n_atoms, double dt, double tau_T, double T_t, double T_eq ){
    alpha_T = 1 + 2*dt / (tau_T) * (T_eq - T_t) / T_t;

    for(int i = 0; i < n_atoms; i++){
        for(int j = 0; j < 3; j++){
            velocities[i][j] *= sqrt(alpha_T);
        }
    }
}

double 
scale_positions(double **positions, int n_atoms, double dt, double tau_P, double P_t, double P_eq, double kappa_T){
    alpha_P = 1 - kappa_T*dt / (tau_P) * (P_eq - P_t);

    for(int i = 0; i < n_atoms; i++){
        for(int j = 0; j < 3; j++){
            positions[i][j] *= cbrt(alpha_P); // cubic root
        }
    }
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
    double P_eq = 1;
    
    double tau_T = dt * 300;
    double tau_P = dt * 300;

    double kappa_T = 1 / (76.E9)

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
    }

    get_forces_AL(forces, positions, N * lattice_param, n_atoms);
    calculate_acceleration(accelerations, forces, n_atoms, masses);
    
    FILE *file = fopen("energies_over_time_dt_0_00001.csv", "w");
    fprintf(file, "time,potential_energy,kinetic_energy,total_energy\n");
    
    printf("Running simulations");
    for (int i = 0; i < total_timesteps; i++) {
        
        velocity_verlet_one_step(accelerations, forces, positions, velocities, dt, lattice_param * N, n_atoms, masses);
        potential_energy = get_energy_AL(positions, N * lattice_param, n_atoms);
        kinetic_energy = calculate_kinetic_energy(velocities, masses, n_atoms);
        total_energy = potential_energy + kinetic_energy;
        
        fprintf(file, "%.4f,%.12f,%.12f,%.12f \n", dt * (1 + i) ,potential_energy, kinetic_energy, total_energy);
        if ((i + 1) % 100 == 0){
            printf("Timestep %d done \n", (1 + i));
        }
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
