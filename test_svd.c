#include "matrix.h"

void test_basic_operations() {
    printf("=== Testing Basic Matrix Operations ===\n");
    
    Matrix* A = matrix_create(3, 3);
    matrix_set(A, 0, 0, 1.0);
    matrix_set(A, 0, 1, 2.0);
    matrix_set(A, 0, 2, 3.0);
    matrix_set(A, 1, 0, 4.0);
    matrix_set(A, 1, 1, 5.0);
    matrix_set(A, 1, 2, 6.0);
    matrix_set(A, 2, 0, 7.0);
    matrix_set(A, 2, 1, 8.0);
    matrix_set(A, 2, 2, 9.0);
    
    printf("Matrix A:\n");
    matrix_print(A);
    
    Matrix* At = matrix_transpose(A);
    printf("A transpose:\n");
    matrix_print(At);
    
    Matrix* AA = matrix_multiply(A, At);
    printf("A * A^T:\n");
    matrix_print(AA);
    
    matrix_destroy(A);
    matrix_destroy(At);
    matrix_destroy(AA);
}

void test_svd_simple() {
    printf("=== Testing SVD on Simple Matrix ===\n");
    
    Matrix* A = matrix_create(2, 2);
    matrix_set(A, 0, 0, 3.0);
    matrix_set(A, 0, 1, 2.0);
    matrix_set(A, 1, 0, 2.0);
    matrix_set(A, 1, 1, 3.0);
    
    printf("Original matrix A:\n");
    matrix_print(A);
    
    SVD* svd = svd_compute(A);
    
    printf("U matrix:\n");
    matrix_print(svd->U);
    
    printf("Singular values (S matrix):\n");
    matrix_print(svd->S);
    
    printf("V matrix:\n");
    matrix_print(svd->V);
    
    printf("Rank: %d\n\n", svd->rank);
    
    Matrix* Vt = matrix_transpose(svd->V);
    Matrix* US = matrix_multiply(svd->U, svd->S);
    Matrix* reconstructed = matrix_multiply(US, Vt);
    
    printf("Reconstructed matrix (U * S * V^T):\n");
    matrix_print(reconstructed);
    
    matrix_destroy(A);
    matrix_destroy(Vt);
    matrix_destroy(US);
    matrix_destroy(reconstructed);
    svd_destroy(svd);
}

void test_svd_rectangular() {
    printf("=== Testing SVD on Rectangular Matrix ===\n");
    
    Matrix* A = matrix_create(3, 2);
    matrix_set(A, 0, 0, 1.0);
    matrix_set(A, 0, 1, 2.0);
    matrix_set(A, 1, 0, 3.0);
    matrix_set(A, 1, 1, 4.0);
    matrix_set(A, 2, 0, 5.0);
    matrix_set(A, 2, 1, 6.0);
    
    printf("Original matrix A (3x2):\n");
    matrix_print(A);
    
    SVD* svd = svd_compute(A);
    
    printf("U matrix:\n");
    matrix_print(svd->U);
    
    printf("Singular values (S matrix):\n");
    matrix_print(svd->S);
    
    printf("V matrix:\n");
    matrix_print(svd->V);
    
    printf("Rank: %d\n\n", svd->rank);
    
    matrix_destroy(A);
    svd_destroy(svd);
}

int main() {
    printf("SVD Implementation Test Program\n");
    printf("===============================\n\n");
    
    test_basic_operations();
    test_svd_simple();
    test_svd_rectangular();
    
    printf("All tests completed!\n");
    return 0;
}