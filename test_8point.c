#include "matrix.h"

void test_essential_matrix_synthetic() {
    printf("=== Testing 8-Point Algorithm with Synthetic Data ===\n");
    
    PointCorrespondences* corr = correspondences_create(8);
    
    correspondences_set(corr, 0, 100, 200, 150, 180);
    correspondences_set(corr, 1, 300, 250, 280, 220);
    correspondences_set(corr, 2, 450, 300, 420, 270);
    correspondences_set(corr, 3, 200, 400, 170, 380);
    correspondences_set(corr, 4, 350, 450, 320, 430);
    correspondences_set(corr, 5, 500, 350, 480, 330);
    correspondences_set(corr, 6, 150, 500, 130, 480);
    correspondences_set(corr, 7, 400, 150, 380, 140);
    
    printf("Point correspondences:\n");
    for (int i = 0; i < corr->count; i++) {
        printf("(%6.1f, %6.1f) <-> (%6.1f, %6.1f)\n",
               corr->points1[i].x, corr->points1[i].y,
               corr->points2[i].x, corr->points2[i].y);
    }
    printf("\n");
    
    Matrix* E = compute_essential_matrix(corr);
    
    printf("Computed Essential Matrix:\n");
    matrix_print(E);
    
    printf("Verification - Essential matrix constraints:\n");
    SVD* svd = svd_compute(E);
    printf("Singular values: %.6f, %.6f, %.6f\n", 
           svd->S->data[0][0], svd->S->data[1][1], svd->S->data[2][2]);
    printf("Rank: %d (should be 2)\n\n", svd->rank);
    
    printf("Epipolar constraint verification (should be ~0):\n");
    for (int i = 0; i < corr->count; i++) {
        double x1 = corr->points1[i].x;
        double y1 = corr->points1[i].y;
        double x2 = corr->points2[i].x;
        double y2 = corr->points2[i].y;
        
        double constraint = x2 * (E->data[0][0] * x1 + E->data[0][1] * y1 + E->data[0][2]) +
                           y2 * (E->data[1][0] * x1 + E->data[1][1] * y1 + E->data[1][2]) +
                                (E->data[2][0] * x1 + E->data[2][1] * y1 + E->data[2][2]);
        
        printf("Point %d: %e\n", i, constraint);
    }
    
    correspondences_destroy(corr);
    matrix_destroy(E);
    svd_destroy(svd);
}

void test_essential_matrix_perfect() {
    printf("=== Testing with Perfect Correspondences ===\n");
    
    PointCorrespondences* corr = correspondences_create(8);
    
    double R[3][3] = {{0.9, -0.1, 0.0}, {0.1, 0.9, 0.0}, {0.0, 0.0, 1.0}};
    double t[3] = {0.1, 0.05, 0.0};
    
    double points_3d[8][3] = {
        {1.0, 0.5, 5.0}, {2.0, 1.0, 4.0}, {-1.0, 1.5, 6.0}, {0.0, -1.0, 3.0},
        {1.5, -0.5, 4.5}, {-0.5, 0.8, 5.5}, {0.8, 1.2, 3.5}, {-1.2, -0.8, 4.8}
    };
    
    for (int i = 0; i < 8; i++) {
        double X = points_3d[i][0];
        double Y = points_3d[i][1];
        double Z = points_3d[i][2];
        
        double x1 = X / Z;
        double y1 = Y / Z;
        
        double X2 = R[0][0] * X + R[0][1] * Y + R[0][2] * Z + t[0];
        double Y2 = R[1][0] * X + R[1][1] * Y + R[1][2] * Z + t[1];
        double Z2 = R[2][0] * X + R[2][1] * Y + R[2][2] * Z + t[2];
        
        double x2 = X2 / Z2;
        double y2 = Y2 / Z2;
        
        correspondences_set(corr, i, x1 * 1000 + 500, y1 * 1000 + 500, 
                                    x2 * 1000 + 500, y2 * 1000 + 500);
    }
    
    printf("Perfect correspondences from known R,t:\n");
    for (int i = 0; i < corr->count; i++) {
        printf("(%6.1f, %6.1f) <-> (%6.1f, %6.1f)\n",
               corr->points1[i].x, corr->points1[i].y,
               corr->points2[i].x, corr->points2[i].y);
    }
    printf("\n");
    
    Matrix* E = compute_essential_matrix(corr);
    
    printf("Computed Essential Matrix:\n");
    matrix_print(E);
    
    printf("Epipolar constraint verification (should be ~0):\n");
    for (int i = 0; i < corr->count; i++) {
        double x1 = corr->points1[i].x;
        double y1 = corr->points1[i].y;
        double x2 = corr->points2[i].x;
        double y2 = corr->points2[i].y;
        
        double constraint = x2 * (E->data[0][0] * x1 + E->data[0][1] * y1 + E->data[0][2]) +
                           y2 * (E->data[1][0] * x1 + E->data[1][1] * y1 + E->data[1][2]) +
                                (E->data[2][0] * x1 + E->data[2][1] * y1 + E->data[2][2]);
        
        printf("Point %d: %e\n", i, constraint);
    }
    
    correspondences_destroy(corr);
    matrix_destroy(E);
}

int main() {
    printf("8-Point Algorithm for Essential Matrix Estimation\n");
    printf("================================================\n\n");
    
    test_essential_matrix_synthetic();
    test_essential_matrix_perfect();
    
    printf("All tests completed!\n");
    return 0;
}