#include "matrix.h"

void test_homography_dlt() {
    printf("=== Testing Homography DLT Estimation ===\n");
    
    PointCorrespondences* corr = correspondences_create(4);
    
    correspondences_set(corr, 0, 0, 0, 100, 50);
    correspondences_set(corr, 1, 200, 0, 300, 80);
    correspondences_set(corr, 2, 200, 150, 320, 200);
    correspondences_set(corr, 3, 0, 150, 80, 180);
    
    printf("Point correspondences:\n");
    for (int i = 0; i < corr->count; i++) {
        printf("(%6.1f, %6.1f) -> (%6.1f, %6.1f)\n",
               corr->points1[i].x, corr->points1[i].y,
               corr->points2[i].x, corr->points2[i].y);
    }
    printf("\n");
    
    Matrix* H = compute_homography_dlt(corr);
    
    printf("Computed Homography Matrix:\n");
    matrix_print(H);
    
    printf("Verification - Forward projection:\n");
    for (int i = 0; i < corr->count; i++) {
        double x1 = corr->points1[i].x;
        double y1 = corr->points1[i].y;
        double x2_expected = corr->points2[i].x;
        double y2_expected = corr->points2[i].y;
        
        double w = H->data[2][0] * x1 + H->data[2][1] * y1 + H->data[2][2];
        double x2_computed = (H->data[0][0] * x1 + H->data[0][1] * y1 + H->data[0][2]) / w;
        double y2_computed = (H->data[1][0] * x1 + H->data[1][1] * y1 + H->data[1][2]) / w;
        
        double error_x = x2_computed - x2_expected;
        double error_y = y2_computed - y2_expected;
        double error = sqrt(error_x * error_x + error_y * error_y);
        
        printf("Point %d: (%6.1f, %6.1f) -> (%6.1f, %6.1f), error: %f\n",
               i, x2_computed, y2_computed, x2_expected, y2_expected, error);
    }
    printf("\n");
    
    correspondences_destroy(corr);
    matrix_destroy(H);
}

void test_homography_decomposition() {
    printf("=== Testing Homography Decomposition ===\n");
    
    Matrix* H = matrix_create(3, 3);
    matrix_set(H, 0, 0, 1.2);
    matrix_set(H, 0, 1, 0.1);
    matrix_set(H, 0, 2, 50.0);
    matrix_set(H, 1, 0, -0.05);
    matrix_set(H, 1, 1, 1.1);
    matrix_set(H, 1, 2, 30.0);
    matrix_set(H, 2, 0, 0.0001);
    matrix_set(H, 2, 1, 0.0002);
    matrix_set(H, 2, 2, 1.0);
    
    printf("Input Homography:\n");
    matrix_print(H);
    
    Matrix* K1 = matrix_create(3, 3);
    matrix_set(K1, 0, 0, 800.0);
    matrix_set(K1, 1, 1, 800.0);
    matrix_set(K1, 0, 2, 320.0);
    matrix_set(K1, 1, 2, 240.0);
    matrix_set(K1, 2, 2, 1.0);
    
    Matrix* K2 = matrix_create(3, 3);
    matrix_set(K2, 0, 0, 800.0);
    matrix_set(K2, 1, 1, 800.0);
    matrix_set(K2, 0, 2, 320.0);
    matrix_set(K2, 1, 2, 240.0);
    matrix_set(K2, 2, 2, 1.0);
    
    printf("Camera intrinsics K:\n");
    matrix_print(K1);
    
    HomographyDecomposition* decomp = decompose_homography(H, K1, K2);
    
    if (decomp) {
        printf("Number of solutions: %d\n\n", decomp->num_solutions);
        
        printf("Solution 1:\n");
        printf("Rotation R1:\n");
        matrix_print(decomp->R1);
        printf("Translation t1:\n");
        matrix_print(decomp->t1);
        printf("Normal n1:\n");
        matrix_print(decomp->n1);
        
        printf("Solution 2:\n");
        printf("Rotation R2:\n");
        matrix_print(decomp->R2);
        printf("Translation t2:\n");
        matrix_print(decomp->t2);
        printf("Normal n2:\n");
        matrix_print(decomp->n2);
        
        double det1 = matrix_determinant_3x3(decomp->R1);
        double det2 = matrix_determinant_3x3(decomp->R2);
        printf("Rotation matrix determinants: %f, %f (should be ~1)\n\n", det1, det2);
        
        homography_decomposition_destroy(decomp);
    } else {
        printf("Decomposition failed!\n");
    }
    
    matrix_destroy(H);
    matrix_destroy(K1);
    matrix_destroy(K2);
}

void test_pseudoinverse() {
    printf("=== Testing Matrix Pseudoinverse ===\n");
    
    Matrix* A = matrix_create(3, 2);
    matrix_set(A, 0, 0, 1.0);
    matrix_set(A, 0, 1, 2.0);
    matrix_set(A, 1, 0, 3.0);
    matrix_set(A, 1, 1, 4.0);
    matrix_set(A, 2, 0, 5.0);
    matrix_set(A, 2, 1, 6.0);
    
    printf("Original matrix A:\n");
    matrix_print(A);
    
    Matrix* pinv = matrix_pseudoinverse(A);
    
    printf("Pseudoinverse A+:\n");
    matrix_print(pinv);
    
    Matrix* At = matrix_transpose(A);
    Matrix* AtA = matrix_multiply(At, A);
    Matrix* should_be_identity = matrix_multiply(pinv, A);
    
    printf("A+ * A (should be identity for full column rank):\n");
    matrix_print(should_be_identity);
    
    matrix_destroy(A);
    matrix_destroy(pinv);
    matrix_destroy(At);
    matrix_destroy(AtA);
    matrix_destroy(should_be_identity);
}

int main() {
    printf("Homography Algorithms Test Program\n");
    printf("==================================\n\n");
    
    test_homography_dlt();
    test_homography_decomposition();
    test_pseudoinverse();
    
    printf("All tests completed!\n");
    return 0;
}