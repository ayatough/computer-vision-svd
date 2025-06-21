#ifndef MATRIX_H
#define MATRIX_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct {
    double **data;
    int rows;
    int cols;
} Matrix;

Matrix* matrix_create(int rows, int cols);
void matrix_destroy(Matrix* m);
void matrix_set(Matrix* m, int row, int col, double value);
double matrix_get(Matrix* m, int row, int col);
void matrix_print(Matrix* m);
Matrix* matrix_transpose(Matrix* m);
Matrix* matrix_multiply(Matrix* a, Matrix* b);
Matrix* matrix_copy(Matrix* m);
void matrix_identity(Matrix* m);

typedef struct {
    Matrix* U;
    Matrix* S;
    Matrix* V;
    int rank;
} SVD;

SVD* svd_compute(Matrix* A);
void svd_destroy(SVD* svd);
void jacobi_eigenvalues(Matrix* A, Matrix* V, double* eigenvals, double tolerance, int max_iterations);

typedef struct {
    double x, y;
} Point2D;

typedef struct {
    Point2D* points1;
    Point2D* points2;
    int count;
} PointCorrespondences;

typedef struct {
    double tx, ty;
    double scale;
} NormalizationParams;

PointCorrespondences* correspondences_create(int count);
void correspondences_destroy(PointCorrespondences* corr);
void correspondences_set(PointCorrespondences* corr, int idx, double x1, double y1, double x2, double y2);
NormalizationParams normalize_points(Point2D* points, int count, Point2D* normalized_points);
void denormalize_matrix(Matrix* E, NormalizationParams norm1, NormalizationParams norm2);
Matrix* build_constraint_matrix(PointCorrespondences* corr);
Matrix* compute_essential_matrix(PointCorrespondences* corr);
void enforce_essential_constraints(Matrix* E);

typedef struct {
    Matrix* R1, *R2;
    Matrix* t1, *t2;
    Matrix* n1, *n2;
    int num_solutions;
} HomographyDecomposition;

Matrix* build_homography_constraint_matrix(PointCorrespondences* corr);
Matrix* compute_homography_dlt(PointCorrespondences* corr);
HomographyDecomposition* decompose_homography(Matrix* H, Matrix* K1, Matrix* K2);
void homography_decomposition_destroy(HomographyDecomposition* decomp);
Matrix* matrix_pseudoinverse(Matrix* A);
double matrix_determinant_3x3(Matrix* m);
Matrix* matrix_inverse_3x3(Matrix* m);

typedef struct {
    Matrix* components;
    Matrix* mean;
    double* eigenvalues;
    double* explained_variance_ratio;
    int n_components;
} PCA;

PCA* pca_fit(Matrix* X, int n_components);
Matrix* pca_transform(PCA* pca, Matrix* X);
Matrix* pca_inverse_transform(PCA* pca, Matrix* X_transformed);
void pca_destroy(PCA* pca);
Matrix* compute_covariance_matrix(Matrix* X);

#endif