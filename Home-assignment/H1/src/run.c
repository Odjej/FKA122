#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "lattice.h"
#include "potential.h"
#include "tools.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

int
run(
    int argc,
    char *argv[]
   )
{
    int N = 4;
    double lattice_param = 5.;
    int n_atoms = 4 * N * N * N;
    
    double **positions = create_2D_array(n_atoms, 3);
    double **forces = create_2D_array(n_atoms, 3);
            
    int total_timesteps = 1000;
    double dt = 0.01;
    init_fcc(positions, N, lattice_param);
    
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    int seed = 42;
    gsl_rng_set(r, seed);

    for(int i = 0; i < n_atoms; i++){
        for(int j = 0; i < 3; j++){
            positions[i][j] += lattice_param * 0.065 * (1- 2 *gsl_rng_uniform(r));
        }
    }

    for(int k = 0; k < total_timesteps; k++){
        get_forces_AL(forces, positions, N * lattice_param, n_atoms);
        
    }

     


    destroy_2D_array(positions, n_atoms);
    return 0;
}
