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
#include "fput.h"
#include "normal_modes.h"

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

#define NPART 32

/* ------------------------------------------------------------------
 * harmonic FPUT chain: all energy initially in mode k = 1
 * ------------------------------------------------------------------ */
static int harmonic_modes(void)
{
    const unsigned int n = NPART;
    const double alpha = 0.0;
    const double E0 = NPART;
    const double dt = 0.1;
    const double t_max = 25000.0;
    const long n_steps = (long) round(t_max / dt);
    const long sample_every = 100;          /* every 10 time units */

    double u[NPART], v[NPART], a[NPART];
    double Q[NPART], P[NPART], energies[NPART];

    /* Q_k = 0, P_k = sqrt(2 E0) delta_k1 */
    for (unsigned int k = 0; k < n; ++k) {
        Q[k] = 0.0;
        P[k] = 0.0;
    }
    P[0] = sqrt(2.0 * E0);

    /* The transform is its own inverse (m = 1) */
    transform_to_normal_modes(Q, u, n);
    transform_to_normal_modes(P, v, n);

    fput_acceleration(a, u, alpha, n);

    FILE *fp = open_output("normal_mode_energies.csv");
    if (fp == NULL)
        return 1;
    fprintf(fp, "time,E1,E2,E3,E4,E5,Etot_modes,Etot_direct\n");

    double max_rel_dev = 0.0;

    for (long step = 0; step <= n_steps; ++step) {
        if (step % sample_every == 0) {
            calculate_normal_mode_energies(energies, u, v, n);

            double e_modes = 0.0;
            for (unsigned int k = 0; k < n; ++k)
                e_modes += energies[k];

            double e_direct = fput_kinetic_energy(v, n)
                            + fput_potential_energy(u, alpha, n);

            double dev = fabs(e_direct - E0) / E0;
            if (dev > max_rel_dev) max_rel_dev = dev;

            fprintf(fp, "%.6f,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e\n",
                    step * dt,
                    energies[0], energies[1], energies[2], energies[3], energies[4],
                    e_modes, e_direct);
        }

        if (step < n_steps)
            fput_verlet_step(a, u, v, alpha, dt, n);
    }
    fclose(fp);

    printf("harmonic modes: max rel. deviation of total energy = %.3e\n", max_rel_dev);
    return 0;
}

/* ------------------------------------------------------------------
 * FPUT chain, general run (tasks 9 and 10)
 *   alpha        : anharmonicity
 *   t_max, dt    : total time and time step
 *   sample_every : steps between samples of the mode energies
 * Writes time, instantaneous E1..E32, running time average A1..A32,
 * and the direct total energy.
 * The time average is the mean over the samples taken so far,
 * which approximates (1/t) * integral of E_k dt'.
 * ------------------------------------------------------------------ */
static int fput_run(double alpha, double t_max, double dt,
                    long sample_every, const char *filename)
{
    const unsigned int n = NPART;
    const double E0 = NPART;
    const long n_steps = (long) round(t_max / dt);

    double u[NPART], v[NPART], a[NPART];
    double Q[NPART], P[NPART];
    double energies[NPART], sum_E[NPART];

    for (unsigned int k = 0; k < n; ++k) {
        Q[k] = 0.0;
        P[k] = 0.0;
        sum_E[k] = 0.0;
    }
    P[0] = sqrt(2.0 * E0);

    transform_to_normal_modes(Q, u, n);
    transform_to_normal_modes(P, v, n);
    fput_acceleration(a, u, alpha, n);

    FILE *fp = open_output(filename);
    if (fp == NULL)
        return 1;

    fprintf(fp, "time");
    for (unsigned int k = 1; k <= n; ++k) fprintf(fp, ",E%u", k);
    for (unsigned int k = 1; k <= n; ++k) fprintf(fp, ",A%u", k);
    fprintf(fp, ",Etot_direct\n");

    long n_samples = 0;

    for (long step = 0; step <= n_steps; ++step) {
        /* step == 0 is written for the instantaneous plot of task 9, but the
         * average starts at the first sample after t = 0 (t_n = n * sample_every * dt) */
        if (step % sample_every == 0) {
            calculate_normal_mode_energies(energies, u, v, n);

            if (step > 0) {
                ++n_samples;
                for (unsigned int k = 0; k < n; ++k)
                    sum_E[k] += energies[k];
            }

            double e_direct = fput_kinetic_energy(v, n)
                            + fput_potential_energy(u, alpha, n);

            fprintf(fp, "%.6f", step * dt);
            for (unsigned int k = 0; k < n; ++k)
                fprintf(fp, ",%.8e", energies[k]);
            for (unsigned int k = 0; k < n; ++k)
                fprintf(fp, ",%.8e", n_samples > 0 ? sum_E[k] / n_samples : energies[k]);
            fprintf(fp, ",%.8e\n", e_direct);
        }

        if (step < n_steps)
            fput_verlet_step(a, u, v, alpha, dt, n);
    }
    fclose(fp);
    printf("fput_run: alpha = %g, t_max = %g, %ld samples written to data/%s\n",
           alpha, t_max, n_samples, filename);
    return 0;
}

/* Task 9: alpha = 0.01 and 0.1, same setup as the harmonic task */
static int anharmonic_short(void)
{
    int status = 0;
    status |= fput_run(0.01, 25000.0, 0.1, 100, "modes_alpha0.01.csv");
    status |= fput_run(0.1,  25000.0, 0.1, 100, "modes_alpha0.1.csv");
    return status;
}

/* Task 10: t_max = 1e6, sample every 1000 steps */
static int anharmonic_long(void)
{
    int status = 0;
    status |= fput_run(0.01, 1.0e6, 0.1, 1000, "long_alpha0.01.csv");
    status |= fput_run(0.1,  1.0e6, 0.1, 1000, "long_alpha0.1.csv");
    return status;
}

/* ------------------------------------------------------------------
 * Dispatcher
 * ------------------------------------------------------------------ */
static void usage(const char *prog)
{
    fprintf(stderr, "Usage: %s [trajectory|scan|spectrum|modes|anharmonic|long|all]\n", prog);
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
    } else if (strcmp(which, "modes") == 0) {
        status = harmonic_modes();
    } else if (strcmp(which, "anharmonic") == 0) {
        status = anharmonic_short();
    } else if (strcmp(which, "long") == 0) {
        status = anharmonic_long();
    } else if (strcmp(which, "all") == 0) {
        status |= trajectory();
        status |= energy_scan();
        status |= spectrum();
        status |= harmonic_modes();
    } else {
        usage(argv[0]);
        return 1;
    }
    return status;
}