#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "unit_conversion.h"
double
amu_to_simulation_units(double mass_amu){
    const double AMU = 1.66053907e-27; // kg
    const double PS = 1e-12; // s
    const double EV = 1.602176634e-19; // J
    const double ANGSTROM = 1e-10; // m 
    double simulation_mass_unit_kg = EV * PS * PS
                                    / (ANGSTROM * ANGSTROM);

    return mass_amu * AMU / simulation_mass_unit_kg;
}

double 
spring_constant_to_simulation_units(double spring_constant){
    const double EV = 1.602176634e-19; // J
    const double ANGSTROM = 1e-10; // m
    double simulation_spring_constant_unit =  EV / (ANGSTROM * ANGSTROM);
    return spring_constant / simulation_spring_constant_unit;
}