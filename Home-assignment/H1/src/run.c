#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "lattice.h"
#include "potential.h"
#include "tools.h"

int
run(
    int argc,
    char *argv[]
   )
{
    int N = 4;
    double lattice_param;
    int n_atoms = 4 * N * N * N;
    
    double **positions = create_2D_array(n_atoms, 3);
    FILE *file = fopen("lattice_potential_energy.csv", "w");
    for(int i = 1; i < 200; i++){
        lattice_param = 3.5 + 0.01 * i;
                
        init_fcc(positions, N, lattice_param);
        
        double E_pot = get_energy_AL(positions, N * lattice_param, n_atoms);
        fprintf(file, "%.6f,%.6f, \n", lattice_param, E_pot);
        
    }
    fclose(file);
    destroy_2D_array(positions, n_atoms);
    return 0;
}
