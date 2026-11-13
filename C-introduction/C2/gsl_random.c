#include "gsl_random.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gsl/gsl_rng.h>

int
write_uniform_sample(size_t count, const char *filename)
{
    double *samples = malloc(count * sizeof(*samples));
    if (samples == NULL) {
        fprintf(stderr, "Could not allocate the uniform sample.\n");
        return -1;
    }

    gsl_rng_env_setup();
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    if (rng == NULL) {
        fprintf(stderr, "Could not allocate the GSL random number generator.\n");
        free(samples);
        return -1;
    }
    gsl_rng_set(rng, (unsigned long)time(NULL));

    for (size_t i = 0; i < count; ++i) {
        samples[i] = gsl_rng_uniform(rng);
    }
    gsl_rng_free(rng);

    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        perror(filename);
        free(samples);
        return -1;
    }

    if (fprintf(file, "uniform\n") < 0) {
        fprintf(stderr, "Could not write the header to %s.\n", filename);
        fclose(file);
        free(samples);
        return -1;
    }
    for (size_t i = 0; i < count; ++i) {
        if (fprintf(file, "%.17g\n", samples[i]) < 0) {
            fprintf(stderr, "Could not write the uniform sample to %s.\n",
                    filename);
            fclose(file);
            free(samples);
            return -1;
        }   
    }

    free(samples);
    if (fclose(file) != 0) {
        perror(filename);
        return -1;
    }
    return 0;
}
