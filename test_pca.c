#include "matrix.h"

void test_pca_2d() {
    printf("=== Testing PCA on 2D Data ===\n");
    
    Matrix* X = matrix_create(6, 2);
    matrix_set(X, 0, 0, 2.5);
    matrix_set(X, 0, 1, 2.4);
    matrix_set(X, 1, 0, 0.5);
    matrix_set(X, 1, 1, 0.7);
    matrix_set(X, 2, 0, 2.2);
    matrix_set(X, 2, 1, 2.9);
    matrix_set(X, 3, 0, 1.9);
    matrix_set(X, 3, 1, 2.2);
    matrix_set(X, 4, 0, 3.1);
    matrix_set(X, 4, 1, 3.0);
    matrix_set(X, 5, 0, 2.3);
    matrix_set(X, 5, 1, 2.7);
    
    printf("Original data:\n");
    matrix_print(X);
    
    PCA* pca = pca_fit(X, 2);
    
    printf("Mean:\n");
    matrix_print(pca->mean);
    
    printf("Principal components:\n");
    matrix_print(pca->components);
    
    printf("Eigenvalues: ");
    for (int i = 0; i < pca->n_components; i++) {
        printf("%.4f ", pca->eigenvalues[i]);
    }
    printf("\n");
    
    printf("Explained variance ratio: ");
    for (int i = 0; i < pca->n_components; i++) {
        printf("%.4f ", pca->explained_variance_ratio[i]);
    }
    printf("\n\n");
    
    Matrix* X_transformed = pca_transform(pca, X);
    printf("Transformed data:\n");
    matrix_print(X_transformed);
    
    Matrix* X_reconstructed = pca_inverse_transform(pca, X_transformed);
    printf("Reconstructed data:\n");
    matrix_print(X_reconstructed);
    
    printf("Reconstruction error:\n");
    double total_error = 0.0;
    for (int i = 0; i < X->rows; i++) {
        for (int j = 0; j < X->cols; j++) {
            double error = X->data[i][j] - X_reconstructed->data[i][j];
            total_error += error * error;
            printf("%8.4f ", error);
        }
        printf("\n");
    }
    printf("Total squared error: %f\n\n", total_error);
    
    matrix_destroy(X);
    matrix_destroy(X_transformed);
    matrix_destroy(X_reconstructed);
    pca_destroy(pca);
}

void test_pca_dimensionality_reduction() {
    printf("=== Testing PCA Dimensionality Reduction ===\n");
    
    Matrix* X = matrix_create(5, 4);
    matrix_set(X, 0, 0, 1.0); matrix_set(X, 0, 1, 2.0); matrix_set(X, 0, 2, 3.0); matrix_set(X, 0, 3, 4.0);
    matrix_set(X, 1, 0, 2.0); matrix_set(X, 1, 1, 4.0); matrix_set(X, 1, 2, 6.0); matrix_set(X, 1, 3, 8.0);
    matrix_set(X, 2, 0, 3.0); matrix_set(X, 2, 1, 6.0); matrix_set(X, 2, 2, 9.0); matrix_set(X, 2, 3, 12.0);
    matrix_set(X, 3, 0, 1.5); matrix_set(X, 3, 1, 3.0); matrix_set(X, 3, 2, 4.5); matrix_set(X, 3, 3, 6.0);
    matrix_set(X, 4, 0, 2.5); matrix_set(X, 4, 1, 5.0); matrix_set(X, 4, 2, 7.5); matrix_set(X, 4, 3, 10.0);
    
    printf("Original 4D data:\n");
    matrix_print(X);
    
    PCA* pca = pca_fit(X, 2);
    
    printf("Reducing to 2 dimensions...\n");
    printf("Explained variance ratio: ");
    for (int i = 0; i < pca->n_components; i++) {
        printf("%.4f ", pca->explained_variance_ratio[i]);
    }
    printf("\n");
    
    double total_explained = 0.0;
    for (int i = 0; i < pca->n_components; i++) {
        total_explained += pca->explained_variance_ratio[i];
    }
    printf("Total explained variance: %.4f (%.1f%%)\n\n", total_explained, total_explained * 100);
    
    Matrix* X_reduced = pca_transform(pca, X);
    printf("Reduced 2D data:\n");
    matrix_print(X_reduced);
    
    Matrix* X_reconstructed = pca_inverse_transform(pca, X_reduced);
    printf("Reconstructed 4D data:\n");
    matrix_print(X_reconstructed);
    
    matrix_destroy(X);
    matrix_destroy(X_reduced);
    matrix_destroy(X_reconstructed);
    pca_destroy(pca);
}

void test_covariance_matrix() {
    printf("=== Testing Covariance Matrix Computation ===\n");
    
    Matrix* X = matrix_create(4, 3);
    matrix_set(X, 0, 0, 1.0); matrix_set(X, 0, 1, 2.0); matrix_set(X, 0, 2, 3.0);
    matrix_set(X, 1, 0, 4.0); matrix_set(X, 1, 1, 5.0); matrix_set(X, 1, 2, 6.0);
    matrix_set(X, 2, 0, 7.0); matrix_set(X, 2, 1, 8.0); matrix_set(X, 2, 2, 9.0);
    matrix_set(X, 3, 0, 2.0); matrix_set(X, 3, 1, 4.0); matrix_set(X, 3, 2, 6.0);
    
    printf("Data matrix:\n");
    matrix_print(X);
    
    Matrix* cov = compute_covariance_matrix(X);
    printf("Covariance matrix:\n");
    matrix_print(cov);
    
    printf("Diagonal elements (variances): ");
    for (int i = 0; i < cov->rows; i++) {
        printf("%.4f ", cov->data[i][i]);
    }
    printf("\n\n");
    
    matrix_destroy(X);
    matrix_destroy(cov);
}

int main() {
    printf("Principal Component Analysis (PCA) Test Program\n");
    printf("===============================================\n\n");
    
    test_covariance_matrix();
    test_pca_2d();
    test_pca_dimensionality_reduction();
    
    printf("All tests completed!\n");
    return 0;
}