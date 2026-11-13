#ifndef LINALG_H
#define LINALG_H


double scalar_product(const double *a, const double *b, int size);
double **allocate_matrix(int rows);
void free_matrix(double **matrix, int rows);
double distance_between_points(const double *point1, const double *point2);
int write_matrix_to_csv(double **matrix, int rows, const char *filename);

#endif
