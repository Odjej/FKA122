#include "linalg.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double
scalar_product(const double *a, const double *b, int size)
{
    double result = 0.0;

    for (int i = 0; i < size; ++i) {
        result += a[i] * b[i];
    }
    return result;
}

double **
allocate_matrix(int rows)
{
    double **matrix = calloc((size_t)rows, sizeof(*matrix));
    if (matrix == NULL) {
        return NULL;
    }

    for (int i = 0; i < rows; ++i) {
        matrix[i] = calloc(3, sizeof(*matrix[i]));
        if (matrix[i] == NULL) {
            free_matrix(matrix, rows);
            return NULL;
        }
    }
    return matrix;
}

void
free_matrix(double **matrix, int rows)
{
    if (matrix == NULL) {
        return;
    }

    for (int i = 0; i < rows; ++i) {
        free(matrix[i]);
    }
    free(matrix);
}

double
distance_between_points(const double *point1, const double *point2)
{
    double squared_distance = 0.0;

    for (int i = 0; i < 3; ++i) {
        double difference = point1[i] - point2[i];
        squared_distance += difference * difference;
    }
    return sqrt(squared_distance);
}

int
write_matrix_to_csv(double **matrix, int rows, const char *filename)
{
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        perror(filename);
        return -1;
    }

    for (int i = 0; i < rows; ++i) {
        if (fprintf(file, "%.6f,%.6f,%.6f\n",
                    matrix[i][0], matrix[i][1], matrix[i][2]) < 0) {
            fprintf(stderr, "Could not write the points matrix to %s.\n",
                    filename);
            fclose(file);
            return -1;
        }
    }

    if (fclose(file) != 0) {
        perror(filename);
        return -1;
    }
    return 0;
}
