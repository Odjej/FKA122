#include <stdio.h>

double
amu_to_simulation_units(double mass_amu)
{
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

/** 
int
main(void)
{
    double carbon_mass_amu = 12.011;
    double oxygen_mass_amu = 15.999;

    double spring_constant = 1600.;

    double carbon_mass = amu_to_simulation_units(carbon_mass_amu);
    double oxygen_mass = amu_to_simulation_units(oxygen_mass_amu);

    double simulation_spring_constant = spring_constant_to_simulation_units(spring_constant);
    
    printf("Carbon mass: %.8f eV*ps^2/Å^2\n", carbon_mass);
    printf("Oxygen mass: %.8f eV*ps^2/Å^2\n", oxygen_mass);
    
    printf("Spring constant: %.8f eV/Å^2\n", simulation_spring_constant);

    return 0;
}
*/
