#include "unit_conversion.h"
#include "tools.h"

int main(int argc, char *argv[]){
    printf("%f\n", amu_to_simulation_units(12.011));
    printf("%f\n", amu_to_simulation_units(15.999));
    printf("%f\n", spring_constant_to_simulation_units(1600.0));
    return 0;
}