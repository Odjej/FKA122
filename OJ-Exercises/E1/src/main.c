#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "unit_conversion.h"
#include "tools.h"
#include "verlet.h"

#define N 3
#define OUTPUT_DIR "data"

/* ------------------------------------------------------------------
 * Opens OUTPUT_DIR/filename for writing, creating the folder if needed.
 * Returns NULL (after printing an error) on failure.
 * ------------------------------------------------------------------ */
static FILE *open_output(const char *filename)
{
    char path[256];

    if (mkdir(OUTPUT_DIR, 0755) != 0 && errno != EEXIST) {
        fprintf(stderr, "Cannot create %s: %s\n", OUTPUT_DIR, strerror(errno));
        return NULL;
    }

    snprintf(path, sizeof path, "%s/%s", OUTPUT_DIR, filename);
    FILE *fp = fopen(path, "w");
    if (fp == NULL)
        fprintf(stderr, "Error opening %s\n", path);
    return fp;
}

/* ------------------------------------------------------------------
 * Shared simulation core
 *
 * Runs n_steps of velocity Verlet from the fixed initial conditions.
 *   fp        : if not NULL, the full trajectory is written as CSV
 *   q1_series : if not NULL, receives q1 at steps 0..n_steps
 * Returns max |E(t) - E(0)| / E(0).
 * ------------------------------------------------------------------ */
static double run(double dt, int n_steps, FILE *fp, double *q1_series)
{
    double masses[N] = {amu_to_simulation_units(15.999),
                        amu_to_simulation_units(12.011),
                        amu_to_simulation_units(15.999)};
    double kappa = spring_constant_to_simulation_units(1600.0);

    double positions[N]  = {0.01, 0.005, -0.005};
    double velocities[N] = {0.0, 0.0, 0.0};
    double accelerations[N];

    calculate_acceleration(accelerations, positions, masses, kappa);

    double e0 = calculate_kinetic_energy(velocities, masses)
              + calculate_potential_energy(positions, kappa);
    double max_dev = 0.0;

    if (fp != NULL)
        fprintf(fp, "time,q1,q2,q3,v1,v2,v3,Ekin,Epot,Etot\n");

    for (int step = 0; step <= n_steps; step++) {
        double ekin = calculate_kinetic_energy(velocities, masses);
        double epot = calculate_potential_energy(positions, kappa);
        double etot = ekin + epot;

        double dev = fabs(etot - e0) / e0;
        if (dev > max_dev) max_dev = dev;

        if (q1_series != NULL) q1_series[step] = positions[0];

        if (fp != NULL) {
            fprintf(fp, "%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10e,%.10e,%.10e\n",
                    step * dt, positions[0], positions[1], positions[2],
                    velocities[0], velocities[1], velocities[2],
                    ekin, epot, etot);
        }

        if (step < n_steps)
            velocity_verlet_one_step(accelerations, positions, velocities,
                                     masses, kappa, dt);
    }
    return max_dev;
}

/* ------------------------------------------------------------------
 * trajectory + energies to CSV 
 * ------------------------------------------------------------------ */
static int trajectory(void)
{
    const double dt = 0.0001;      /* ps */
    const int n_steps = 10000;

    FILE *fp = open_output("verlet_output.csv");
    if (fp == NULL)
        return 1;

    double dev = run(dt, n_steps, fp, NULL);
    fclose(fp);

    printf("trajectory: dt = %g ps, max rel. energy deviation = %.3e\n", dt, dev);
    return 0;
}

/* ------------------------------------------------------------------
 *  energy conservation vs time step
 * ------------------------------------------------------------------ */
static int energy_scan(void)
{
    double dts[] = {0.002, 0.001, 0.0005, 0.0002, 0.0001, 0.00005};
    int n_dts = sizeof(dts) / sizeof(dts[0]);

    FILE *fe = open_output("energy_vs_dt.csv");
    if (fe == NULL)
        return 1;

    fprintf(fe, "dt_ps,max_rel_dev\n");
    for (int k = 0; k < n_dts; k++) {
        int steps = (int)round(1.0 / dts[k]);
        double dev = run(dts[k], steps, NULL, NULL);
        printf("scan: dt = %.5f ps  max rel. deviation = %.3e\n", dts[k], dev);
        fprintf(fe, "%.6f,%.6e\n", dts[k], dev);
    }
    fclose(fe);
    return 0;
}

/* ------------------------------------------------------------------
 * power spectrum
 * ------------------------------------------------------------------ */
static int spectrum(void)
{
    const double dt = 0.0001;      /* ps, chosen from energy scan */
    const int n = 100000;

    double *q1    = malloc(sizeof(double) * n);
    double *power = malloc(sizeof(double) * n);
    double *omega = malloc(sizeof(double) * n);
    if (q1 == NULL || power == NULL || omega == NULL) {
        fprintf(stderr, "Out of memory\n");
        free(q1); free(power); free(omega);
        return 1;
    }

    run(dt, n - 1, NULL, q1);

    double mean = average(q1, n);
    addition_with_constant(q1, q1, -mean, n);   /* optional: drop the DC mode */
    hann_window(q1, q1, n);                     /* optional: suppress leakage */

    powerspectrum(power, q1, n, dt);
    multiplication_with_constant(power, power, 1.0 / hann_power_correction(n), n);
    fft_freq(omega, n, dt);

    FILE *fs = open_output("spectrum.csv");
    if (fs == NULL) {
        free(q1); free(power); free(omega);
        return 1;
    }
    fprintf(fs, "omega_rad_per_ps,freq_THz,power\n");
    double pmax = 0.0;
    for (int i = 0; i <= n / 2; i++) {
        fprintf(fs, "%.8e,%.8e,%.8e\n", omega[i], omega[i] / (2.0 * M_PI), power[i]);
        if (i > 0 && power[i] > pmax) pmax = power[i];
    }
    fclose(fs);

    printf("spectrum peaks:\n");
    if (power[0] > power[1])
        printf("  omega = %.2f rad/ps   f = %.3f THz\n", omega[0], omega[0] / (2.0 * M_PI));
    for (int i = 1; i < n / 2; i++) {
        if (power[i] > power[i - 1] && power[i] > power[i + 1]
            && power[i] > 1e-3 * pmax)
            printf("  omega = %.2f rad/ps   f = %.3f THz\n",
                omega[i], omega[i] / (2.0 * M_PI));
    }

    free(q1); free(power); free(omega);
    return 0;

    // Reference values from experimental values for the CO2 bending mode:
    // 1351cm^-1 = 40.5 THz = 254 rad/ps
    // 2349cm^-1 = 70.4 THz = 442 rad/ps
    // 667cm^-1 = 20.0 THz = 126 rad/ps, but this is bend not stetch, so we don't see it here.
}

/* ------------------------------------------------------------------
 * Dispatcher
 * ------------------------------------------------------------------ */
static void usage(const char *prog)
{
    fprintf(stderr, "Usage: %s [trajectory|scan|spectrum|all]\n", prog);
}

int main(int argc, char *argv[])
{
    const char *which = (argc > 1) ? argv[1] : "all";
    int status = 0;

    if (strcmp(which, "trajectory") == 0) {
        status = trajectory();
    } else if (strcmp(which, "scan") == 0) {
        status = energy_scan();    
    } else if (strcmp(which, "spectrum") == 0) {
        status = spectrum();
    } else if (strcmp(which, "all") == 0) {
        status |= trajectory();
        status |= energy_scan();
        status |= spectrum();
    } else {
        usage(argv[0]);
        return 1;
    }
    return status;
}