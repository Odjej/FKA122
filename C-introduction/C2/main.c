#include "gsl_random.h"
#include "linalg.h"

#include <stdio.h>
#include <stdlib.h>

#define RANDOM_SAMPLE_SIZE 200000

int
main(void)
{
    int size;
    double *a;
    double *b;
    double **points;

    printf("Enter the number of points (at least 2): ");
    if (scanf("%d", &size) != 1 || size < 2) {
        fprintf(stderr, "Please enter an integer of 2 or greater.\n");
        return EXIT_FAILURE;
    }

    a = calloc((size_t)size, sizeof(*a));
    b = calloc((size_t)size, sizeof(*b));
    if (a == NULL || b == NULL) {
        fprintf(stderr, "Could not allocate the vectors.\n");
        free(a);
        free(b);
        return EXIT_FAILURE;
    }

    for (int i = 0; i < size; ++i) {
        a[i] = (double)(i + 1);
        b[i] = (double)(i + 2);
    }
    printf("Scalar product: %.2f\n", scalar_product(a, b, size));
    free(a);
    free(b);

    points = allocate_matrix(size);
    if (points == NULL) {
        fprintf(stderr, "Could not allocate the points matrix.\n");
        return EXIT_FAILURE;
    }

    for (int i = 0; i < size; ++i) {
        points[i][0] = (double)i;
        points[i][1] = (double)(i * 2);
        points[i][2] = (double)(i * 3);
    }

    printf("Distance between the first two points: %.2f\n",
           distance_between_points(points[0], points[1]));
    if (write_matrix_to_csv(points, size, "points.csv") != 0) {
        free_matrix(points, size);
        return EXIT_FAILURE;
    }
    printf("Wrote %d points to points.csv\n", size);

    if (write_uniform_sample(RANDOM_SAMPLE_SIZE, "uniform_samples.csv") != 0) {
        free_matrix(points, size);
        return EXIT_FAILURE;
    }
    printf("Wrote %d uniform samples to uniform_samples.csv\n",
           RANDOM_SAMPLE_SIZE);

    free_matrix(points, size);
    return EXIT_SUCCESS;
}
