#pragma once
#include <stdio.h>

// Function declarations for unit conversion

// Converts mass from atomic mass units (amu) to simulation units (eV*ps^2/Å^2)
double
amu_to_simulation_units(double mass_amu);

// Converts spring constant from N/m to simulation units (eV/Å^2)
double
spring_constant_to_simulation_units(double spring_constant);