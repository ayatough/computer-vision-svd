#include "matrix.h"

Matrix* matrix_create(int rows, int cols) {
    Matrix* m = (Matrix*)malloc(sizeof(Matrix));
    if (!m) return NULL;
    
    m->rows = rows;
    m->cols = cols;
    m->data = (double**)malloc(rows * sizeof(double*));
    if (!m->data) {
        free(m);
        return NULL;
    }
    
    for (int i = 0; i < rows; i++) {
        m->data[i] = (double*)calloc(cols, sizeof(double));
        if (!m->data[i]) {
            for (int j = 0; j < i; j++) {
                free(m->data[j]);
            }
            free(m->data);
            free(m);
            return NULL;
        }
    }
    
    return m;
}

void matrix_destroy(Matrix* m) {
    if (!m) return;
    
    for (int i = 0; i < m->rows; i++) {
        free(m->data[i]);
    }
    free(m->data);
    free(m);
}

void matrix_set(Matrix* m, int row, int col, double value) {
    if (m && row >= 0 && row < m->rows && col >= 0 && col < m->cols) {
        m->data[row][col] = value;
    }
}

double matrix_get(Matrix* m, int row, int col) {
    if (m && row >= 0 && row < m->rows && col >= 0 && col < m->cols) {
        return m->data[row][col];
    }
    return 0.0;
}

void matrix_print(Matrix* m) {
    if (!m) return;
    
    for (int i = 0; i < m->rows; i++) {
        for (int j = 0; j < m->cols; j++) {
            printf("%8.4f ", m->data[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

Matrix* matrix_copy(Matrix* m) {
    if (!m) return NULL;
    
    Matrix* copy = matrix_create(m->rows, m->cols);
    if (!copy) return NULL;
    
    for (int i = 0; i < m->rows; i++) {
        for (int j = 0; j < m->cols; j++) {
            copy->data[i][j] = m->data[i][j];
        }
    }
    
    return copy;
}

void matrix_identity(Matrix* m) {
    if (!m) return;
    
    for (int i = 0; i < m->rows; i++) {
        for (int j = 0; j < m->cols; j++) {
            m->data[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }
}

Matrix* matrix_transpose(Matrix* m) {
    if (!m) return NULL;
    
    Matrix* transpose = matrix_create(m->cols, m->rows);
    if (!transpose) return NULL;
    
    for (int i = 0; i < m->rows; i++) {
        for (int j = 0; j < m->cols; j++) {
            transpose->data[j][i] = m->data[i][j];
        }
    }
    
    return transpose;
}

Matrix* matrix_multiply(Matrix* a, Matrix* b) {
    if (!a || !b || a->cols != b->rows) return NULL;
    
    Matrix* result = matrix_create(a->rows, b->cols);
    if (!result) return NULL;
    
    for (int i = 0; i < a->rows; i++) {
        for (int j = 0; j < b->cols; j++) {
            double sum = 0.0;
            for (int k = 0; k < a->cols; k++) {
                sum += a->data[i][k] * b->data[k][j];
            }
            result->data[i][j] = sum;
        }
    }
    
    return result;
}

void jacobi_eigenvalues(Matrix* A, Matrix* V, double* eigenvals, double tolerance, int max_iterations) {
    if (!A || !V || !eigenvals) return;
    
    int n = A->rows;
    Matrix* B = matrix_copy(A);
    matrix_identity(V);
    
    for (int iter = 0; iter < max_iterations; iter++) {
        double max_off_diag = 0.0;
        int p = 0, q = 1;
        
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                double val = fabs(B->data[i][j]);
                if (val > max_off_diag) {
                    max_off_diag = val;
                    p = i;
                    q = j;
                }
            }
        }
        
        if (max_off_diag < tolerance) break;
        
        double app = B->data[p][p];
        double aqq = B->data[q][q];
        double apq = B->data[p][q];
        
        double theta = (aqq - app) / (2.0 * apq);
        double t = (theta >= 0) ? 1.0 / (theta + sqrt(1.0 + theta * theta)) : 1.0 / (theta - sqrt(1.0 + theta * theta));
        double c = 1.0 / sqrt(1.0 + t * t);
        double s = t * c;
        
        for (int i = 0; i < n; i++) {
            if (i != p && i != q) {
                double bip = B->data[i][p];
                double biq = B->data[i][q];
                B->data[i][p] = c * bip - s * biq;
                B->data[p][i] = B->data[i][p];
                B->data[i][q] = s * bip + c * biq;
                B->data[q][i] = B->data[i][q];
            }
        }
        
        B->data[p][p] = c * c * app + s * s * aqq - 2.0 * c * s * apq;
        B->data[q][q] = s * s * app + c * c * aqq + 2.0 * c * s * apq;
        B->data[p][q] = 0.0;
        B->data[q][p] = 0.0;
        
        for (int i = 0; i < n; i++) {
            double vip = V->data[i][p];
            double viq = V->data[i][q];
            V->data[i][p] = c * vip - s * viq;
            V->data[i][q] = s * vip + c * viq;
        }
    }
    
    for (int i = 0; i < n; i++) {
        eigenvals[i] = B->data[i][i];
    }
    
    for (int i = 0; i < n - 1; i++) {
        for (int j = i + 1; j < n; j++) {
            if (eigenvals[i] < eigenvals[j]) {
                double temp = eigenvals[i];
                eigenvals[i] = eigenvals[j];
                eigenvals[j] = temp;
                
                for (int k = 0; k < n; k++) {
                    temp = V->data[k][i];
                    V->data[k][i] = V->data[k][j];
                    V->data[k][j] = temp;
                }
            }
        }
    }
    
    matrix_destroy(B);
}

SVD* svd_compute(Matrix* A) {
    if (!A) return NULL;
    
    SVD* svd = (SVD*)malloc(sizeof(SVD));
    if (!svd) return NULL;
    
    int m = A->rows;
    int n = A->cols;
    int min_dim = (m < n) ? m : n;
    
    Matrix* At = matrix_transpose(A);
    Matrix* AtA = matrix_multiply(At, A);
    Matrix* AAt = matrix_multiply(A, At);
    
    Matrix* V = matrix_create(n, n);
    Matrix* U = matrix_create(m, m);
    double* sigma_squared = (double*)malloc(min_dim * sizeof(double));
    double* sigma_squared_u = (double*)malloc(m * sizeof(double));
    
    jacobi_eigenvalues(AtA, V, sigma_squared, 1e-12, 1000);
    jacobi_eigenvalues(AAt, U, sigma_squared_u, 1e-12, 1000);
    
    svd->S = matrix_create(m, n);
    svd->U = matrix_create(m, m);
    svd->V = matrix_create(n, n);
    
    int rank = 0;
    for (int i = 0; i < min_dim; i++) {
        if (sigma_squared[i] > 1e-12) {
            svd->S->data[i][i] = sqrt(sigma_squared[i]);
            rank++;
        }
    }
    svd->rank = rank;
    
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            svd->U->data[i][j] = U->data[i][j];
        }
    }
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            svd->V->data[i][j] = V->data[i][j];
        }
    }
    
    for (int i = 0; i < rank; i++) {
        Matrix* Av = matrix_create(m, 1);
        for (int j = 0; j < m; j++) {
            double sum = 0.0;
            for (int k = 0; k < n; k++) {
                sum += A->data[j][k] * svd->V->data[k][i];
            }
            Av->data[j][0] = sum;
        }
        
        double norm = 0.0;
        for (int j = 0; j < m; j++) {
            norm += Av->data[j][0] * Av->data[j][0];
        }
        norm = sqrt(norm);
        
        if (norm > 1e-12) {
            for (int j = 0; j < m; j++) {
                svd->U->data[j][i] = Av->data[j][0] / norm;
            }
        }
        
        matrix_destroy(Av);
    }
    
    matrix_destroy(At);
    matrix_destroy(AtA);
    matrix_destroy(AAt);
    matrix_destroy(V);
    matrix_destroy(U);
    free(sigma_squared);
    free(sigma_squared_u);
    
    return svd;
}

void svd_destroy(SVD* svd) {
    if (!svd) return;
    
    matrix_destroy(svd->U);
    matrix_destroy(svd->S);
    matrix_destroy(svd->V);
    free(svd);
}

PointCorrespondences* correspondences_create(int count) {
    PointCorrespondences* corr = (PointCorrespondences*)malloc(sizeof(PointCorrespondences));
    if (!corr) return NULL;
    
    corr->count = count;
    corr->points1 = (Point2D*)malloc(count * sizeof(Point2D));
    corr->points2 = (Point2D*)malloc(count * sizeof(Point2D));
    
    if (!corr->points1 || !corr->points2) {
        free(corr->points1);
        free(corr->points2);
        free(corr);
        return NULL;
    }
    
    return corr;
}

void correspondences_destroy(PointCorrespondences* corr) {
    if (!corr) return;
    
    free(corr->points1);
    free(corr->points2);
    free(corr);
}

void correspondences_set(PointCorrespondences* corr, int idx, double x1, double y1, double x2, double y2) {
    if (!corr || idx < 0 || idx >= corr->count) return;
    
    corr->points1[idx].x = x1;
    corr->points1[idx].y = y1;
    corr->points2[idx].x = x2;
    corr->points2[idx].y = y2;
}

NormalizationParams normalize_points(Point2D* points, int count, Point2D* normalized_points) {
    NormalizationParams params;
    
    double cx = 0.0, cy = 0.0;
    for (int i = 0; i < count; i++) {
        cx += points[i].x;
        cy += points[i].y;
    }
    cx /= count;
    cy /= count;
    
    double avg_dist = 0.0;
    for (int i = 0; i < count; i++) {
        double dx = points[i].x - cx;
        double dy = points[i].y - cy;
        avg_dist += sqrt(dx * dx + dy * dy);
    }
    avg_dist /= count;
    
    double scale = sqrt(2.0) / avg_dist;
    
    for (int i = 0; i < count; i++) {
        normalized_points[i].x = scale * (points[i].x - cx);
        normalized_points[i].y = scale * (points[i].y - cy);
    }
    
    params.tx = cx;
    params.ty = cy;
    params.scale = scale;
    
    return params;
}

void denormalize_matrix(Matrix* E, NormalizationParams norm1, NormalizationParams norm2) {
    if (!E) return;
    
    Matrix* T1 = matrix_create(3, 3);
    Matrix* T2 = matrix_create(3, 3);
    
    matrix_set(T1, 0, 0, norm1.scale);
    matrix_set(T1, 1, 1, norm1.scale);
    matrix_set(T1, 2, 2, 1.0);
    matrix_set(T1, 0, 2, -norm1.scale * norm1.tx);
    matrix_set(T1, 1, 2, -norm1.scale * norm1.ty);
    
    matrix_set(T2, 0, 0, norm2.scale);
    matrix_set(T2, 1, 1, norm2.scale);
    matrix_set(T2, 2, 2, 1.0);
    matrix_set(T2, 0, 2, -norm2.scale * norm2.tx);
    matrix_set(T2, 1, 2, -norm2.scale * norm2.ty);
    
    Matrix* T2t = matrix_transpose(T2);
    Matrix* temp = matrix_multiply(T2t, E);
    Matrix* result = matrix_multiply(temp, T1);
    
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            E->data[i][j] = result->data[i][j];
        }
    }
    
    matrix_destroy(T1);
    matrix_destroy(T2);
    matrix_destroy(T2t);
    matrix_destroy(temp);
    matrix_destroy(result);
}

Matrix* build_constraint_matrix(PointCorrespondences* corr) {
    if (!corr || corr->count < 8) return NULL;
    
    Matrix* A = matrix_create(corr->count, 9);
    
    for (int i = 0; i < corr->count; i++) {
        double x1 = corr->points1[i].x;
        double y1 = corr->points1[i].y;
        double x2 = corr->points2[i].x;
        double y2 = corr->points2[i].y;
        
        matrix_set(A, i, 0, x1 * x2);
        matrix_set(A, i, 1, x1 * y2);
        matrix_set(A, i, 2, x1);
        matrix_set(A, i, 3, y1 * x2);
        matrix_set(A, i, 4, y1 * y2);
        matrix_set(A, i, 5, y1);
        matrix_set(A, i, 6, x2);
        matrix_set(A, i, 7, y2);
        matrix_set(A, i, 8, 1.0);
    }
    
    return A;
}

void enforce_essential_constraints(Matrix* E) {
    if (!E || E->rows != 3 || E->cols != 3) return;
    
    SVD* svd = svd_compute(E);
    
    double s1 = svd->S->data[0][0];
    double s2 = svd->S->data[1][1];
    double avg_singular = (s1 + s2) / 2.0;
    
    svd->S->data[0][0] = avg_singular;
    svd->S->data[1][1] = avg_singular;
    svd->S->data[2][2] = 0.0;
    
    Matrix* Vt = matrix_transpose(svd->V);
    Matrix* US = matrix_multiply(svd->U, svd->S);
    Matrix* corrected = matrix_multiply(US, Vt);
    
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            E->data[i][j] = corrected->data[i][j];
        }
    }
    
    matrix_destroy(Vt);
    matrix_destroy(US);
    matrix_destroy(corrected);
    svd_destroy(svd);
}

Matrix* compute_essential_matrix(PointCorrespondences* corr) {
    if (!corr || corr->count < 8) return NULL;
    
    Point2D* norm_points1 = (Point2D*)malloc(corr->count * sizeof(Point2D));
    Point2D* norm_points2 = (Point2D*)malloc(corr->count * sizeof(Point2D));
    
    NormalizationParams norm1 = normalize_points(corr->points1, corr->count, norm_points1);
    NormalizationParams norm2 = normalize_points(corr->points2, corr->count, norm_points2);
    
    PointCorrespondences* norm_corr = correspondences_create(corr->count);
    for (int i = 0; i < corr->count; i++) {
        correspondences_set(norm_corr, i, 
                          norm_points1[i].x, norm_points1[i].y,
                          norm_points2[i].x, norm_points2[i].y);
    }
    
    Matrix* A = build_constraint_matrix(norm_corr);
    SVD* svd = svd_compute(A);
    
    Matrix* E = matrix_create(3, 3);
    for (int i = 0; i < 9; i++) {
        int row = i / 3;
        int col = i % 3;
        E->data[row][col] = svd->V->data[i][svd->V->cols - 1];
    }
    
    enforce_essential_constraints(E);
    denormalize_matrix(E, norm1, norm2);
    
    free(norm_points1);
    free(norm_points2);
    correspondences_destroy(norm_corr);
    matrix_destroy(A);
    svd_destroy(svd);
    
    return E;
}

Matrix* matrix_pseudoinverse(Matrix* A) {
    if (!A) return NULL;
    
    SVD* svd = svd_compute(A);
    if (!svd) return NULL;
    
    Matrix* S_inv = matrix_create(svd->S->cols, svd->S->rows);
    
    for (int i = 0; i < svd->S->rows && i < svd->S->cols; i++) {
        double s = svd->S->data[i][i];
        if (s > 1e-12) {
            S_inv->data[i][i] = 1.0 / s;
        }
    }
    
    Matrix* Ut = matrix_transpose(svd->U);
    Matrix* temp = matrix_multiply(svd->V, S_inv);
    Matrix* pinv = matrix_multiply(temp, Ut);
    
    matrix_destroy(S_inv);
    matrix_destroy(Ut);
    matrix_destroy(temp);
    svd_destroy(svd);
    
    return pinv;
}

double matrix_determinant_3x3(Matrix* m) {
    if (!m || m->rows != 3 || m->cols != 3) return 0.0;
    
    double det = m->data[0][0] * (m->data[1][1] * m->data[2][2] - m->data[1][2] * m->data[2][1])
               - m->data[0][1] * (m->data[1][0] * m->data[2][2] - m->data[1][2] * m->data[2][0])
               + m->data[0][2] * (m->data[1][0] * m->data[2][1] - m->data[1][1] * m->data[2][0]);
    
    return det;
}

Matrix* matrix_inverse_3x3(Matrix* m) {
    if (!m || m->rows != 3 || m->cols != 3) return NULL;
    
    double det = matrix_determinant_3x3(m);
    if (fabs(det) < 1e-12) return NULL;
    
    Matrix* inv = matrix_create(3, 3);
    
    inv->data[0][0] = (m->data[1][1] * m->data[2][2] - m->data[1][2] * m->data[2][1]) / det;
    inv->data[0][1] = (m->data[0][2] * m->data[2][1] - m->data[0][1] * m->data[2][2]) / det;
    inv->data[0][2] = (m->data[0][1] * m->data[1][2] - m->data[0][2] * m->data[1][1]) / det;
    
    inv->data[1][0] = (m->data[1][2] * m->data[2][0] - m->data[1][0] * m->data[2][2]) / det;
    inv->data[1][1] = (m->data[0][0] * m->data[2][2] - m->data[0][2] * m->data[2][0]) / det;
    inv->data[1][2] = (m->data[0][2] * m->data[1][0] - m->data[0][0] * m->data[1][2]) / det;
    
    inv->data[2][0] = (m->data[1][0] * m->data[2][1] - m->data[1][1] * m->data[2][0]) / det;
    inv->data[2][1] = (m->data[0][1] * m->data[2][0] - m->data[0][0] * m->data[2][1]) / det;
    inv->data[2][2] = (m->data[0][0] * m->data[1][1] - m->data[0][1] * m->data[1][0]) / det;
    
    return inv;
}

Matrix* build_homography_constraint_matrix(PointCorrespondences* corr) {
    if (!corr || corr->count < 4) return NULL;
    
    Matrix* A = matrix_create(2 * corr->count, 9);
    
    for (int i = 0; i < corr->count; i++) {
        double x1 = corr->points1[i].x;
        double y1 = corr->points1[i].y;
        double x2 = corr->points2[i].x;
        double y2 = corr->points2[i].y;
        
        int row1 = 2 * i;
        int row2 = 2 * i + 1;
        
        matrix_set(A, row1, 0, -x1);
        matrix_set(A, row1, 1, -y1);
        matrix_set(A, row1, 2, -1.0);
        matrix_set(A, row1, 3, 0.0);
        matrix_set(A, row1, 4, 0.0);
        matrix_set(A, row1, 5, 0.0);
        matrix_set(A, row1, 6, x2 * x1);
        matrix_set(A, row1, 7, x2 * y1);
        matrix_set(A, row1, 8, x2);
        
        matrix_set(A, row2, 0, 0.0);
        matrix_set(A, row2, 1, 0.0);
        matrix_set(A, row2, 2, 0.0);
        matrix_set(A, row2, 3, -x1);
        matrix_set(A, row2, 4, -y1);
        matrix_set(A, row2, 5, -1.0);
        matrix_set(A, row2, 6, y2 * x1);
        matrix_set(A, row2, 7, y2 * y1);
        matrix_set(A, row2, 8, y2);
    }
    
    return A;
}

Matrix* compute_homography_dlt(PointCorrespondences* corr) {
    if (!corr || corr->count < 4) return NULL;
    
    Point2D* norm_points1 = (Point2D*)malloc(corr->count * sizeof(Point2D));
    Point2D* norm_points2 = (Point2D*)malloc(corr->count * sizeof(Point2D));
    
    NormalizationParams norm1 = normalize_points(corr->points1, corr->count, norm_points1);
    NormalizationParams norm2 = normalize_points(corr->points2, corr->count, norm_points2);
    
    PointCorrespondences* norm_corr = correspondences_create(corr->count);
    for (int i = 0; i < corr->count; i++) {
        correspondences_set(norm_corr, i, 
                          norm_points1[i].x, norm_points1[i].y,
                          norm_points2[i].x, norm_points2[i].y);
    }
    
    Matrix* A = build_homography_constraint_matrix(norm_corr);
    SVD* svd = svd_compute(A);
    
    Matrix* H = matrix_create(3, 3);
    for (int i = 0; i < 9; i++) {
        int row = i / 3;
        int col = i % 3;
        H->data[row][col] = svd->V->data[i][svd->V->cols - 1];
    }
    
    Matrix* T1 = matrix_create(3, 3);
    Matrix* T2 = matrix_create(3, 3);
    
    matrix_set(T1, 0, 0, norm1.scale);
    matrix_set(T1, 1, 1, norm1.scale);
    matrix_set(T1, 2, 2, 1.0);
    matrix_set(T1, 0, 2, -norm1.scale * norm1.tx);
    matrix_set(T1, 1, 2, -norm1.scale * norm1.ty);
    
    matrix_set(T2, 0, 0, norm2.scale);
    matrix_set(T2, 1, 1, norm2.scale);
    matrix_set(T2, 2, 2, 1.0);
    matrix_set(T2, 0, 2, -norm2.scale * norm2.tx);
    matrix_set(T2, 1, 2, -norm2.scale * norm2.ty);
    
    Matrix* T2_inv = matrix_inverse_3x3(T2);
    Matrix* temp = matrix_multiply(T2_inv, H);
    Matrix* final_H = matrix_multiply(temp, T1);
    
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            H->data[i][j] = final_H->data[i][j];
        }
    }
    
    free(norm_points1);
    free(norm_points2);
    correspondences_destroy(norm_corr);
    matrix_destroy(A);
    matrix_destroy(T1);
    matrix_destroy(T2);
    matrix_destroy(T2_inv);
    matrix_destroy(temp);
    matrix_destroy(final_H);
    svd_destroy(svd);
    
    return H;
}

HomographyDecomposition* decompose_homography(Matrix* H, Matrix* K1, Matrix* K2) {
    if (!H || H->rows != 3 || H->cols != 3) return NULL;
    
    HomographyDecomposition* decomp = (HomographyDecomposition*)malloc(sizeof(HomographyDecomposition));
    if (!decomp) return NULL;
    
    Matrix* K1_inv = matrix_inverse_3x3(K1);
    Matrix* K2_inv = matrix_inverse_3x3(K2);
    
    if (!K1_inv || !K2_inv) {
        matrix_destroy(K1_inv);
        matrix_destroy(K2_inv);
        free(decomp);
        return NULL;
    }
    
    Matrix* temp = matrix_multiply(K2_inv, H);
    Matrix* A = matrix_multiply(temp, K1);
    
    SVD* svd = svd_compute(A);
    
    double sigma1 = svd->S->data[0][0];
    double sigma2 = svd->S->data[1][1];
    double sigma3 = svd->S->data[2][2];
    
    decomp->R1 = matrix_create(3, 3);
    decomp->R2 = matrix_create(3, 3);
    decomp->t1 = matrix_create(3, 1);
    decomp->t2 = matrix_create(3, 1);
    decomp->n1 = matrix_create(3, 1);
    decomp->n2 = matrix_create(3, 1);
    
    Matrix* U = svd->U;
    Matrix* V = svd->V;
    Matrix* Vt = matrix_transpose(V);
    
    double d1 = sigma1;
    double d2 = sigma3;
    double d3 = sigma2;
    
    if (d1 < d3) {
        double temp_d = d1; d1 = d3; d3 = temp_d;
    }
    
    double cos_theta = (d1*d1 + d3*d3 - d2*d2) / (2*d1*d3);
    if (cos_theta > 1.0) cos_theta = 1.0;
    if (cos_theta < -1.0) cos_theta = -1.0;
    double sin_theta = sqrt(1.0 - cos_theta*cos_theta);
    
    double n1_vals[4][3] = {
        {sin_theta, 0, cos_theta},
        {-sin_theta, 0, -cos_theta},
        {sin_theta, 0, cos_theta},
        {-sin_theta, 0, -cos_theta}
    };
    
    for (int sol = 0; sol < 2; sol++) {
        Matrix* R_temp = matrix_create(3, 3);
        Matrix* t_temp = matrix_create(3, 1);
        Matrix* n_temp = matrix_create(3, 1);
        
        for (int i = 0; i < 3; i++) {
            n_temp->data[i][0] = n1_vals[sol][i];
        }
        
        Matrix* W = matrix_create(3, 3);
        matrix_set(W, 0, 0, cos_theta);
        matrix_set(W, 0, 2, sin_theta * ((sol == 0) ? 1 : -1));
        matrix_set(W, 1, 1, 1.0);
        matrix_set(W, 2, 0, -sin_theta * ((sol == 0) ? 1 : -1));
        matrix_set(W, 2, 2, cos_theta);
        
        Matrix* UW = matrix_multiply(U, W);
        Matrix* R = matrix_multiply(UW, Vt);
        
        double t_scale = (d1 - d3) / d2;
        for (int i = 0; i < 3; i++) {
            t_temp->data[i][0] = t_scale * U->data[i][1];
        }
        
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                R_temp->data[i][j] = R->data[i][j];
            }
        }
        
        if (sol == 0) {
            decomp->R1 = R_temp;
            decomp->t1 = t_temp;
            decomp->n1 = n_temp;
        } else {
            decomp->R2 = R_temp;
            decomp->t2 = t_temp;
            decomp->n2 = n_temp;
        }
        
        matrix_destroy(W);
        matrix_destroy(UW);
        matrix_destroy(R);
    }
    
    decomp->num_solutions = 2;
    
    matrix_destroy(K1_inv);
    matrix_destroy(K2_inv);
    matrix_destroy(temp);
    matrix_destroy(A);
    matrix_destroy(Vt);
    svd_destroy(svd);
    
    return decomp;
}

void homography_decomposition_destroy(HomographyDecomposition* decomp) {
    if (!decomp) return;
    
    matrix_destroy(decomp->R1);
    matrix_destroy(decomp->R2);
    matrix_destroy(decomp->t1);
    matrix_destroy(decomp->t2);
    matrix_destroy(decomp->n1);
    matrix_destroy(decomp->n2);
    free(decomp);
}

Matrix* compute_covariance_matrix(Matrix* X) {
    if (!X) return NULL;
    
    int n_samples = X->rows;
    int n_features = X->cols;
    
    Matrix* mean = matrix_create(1, n_features);
    for (int j = 0; j < n_features; j++) {
        double sum = 0.0;
        for (int i = 0; i < n_samples; i++) {
            sum += X->data[i][j];
        }
        mean->data[0][j] = sum / n_samples;
    }
    
    Matrix* X_centered = matrix_create(n_samples, n_features);
    for (int i = 0; i < n_samples; i++) {
        for (int j = 0; j < n_features; j++) {
            X_centered->data[i][j] = X->data[i][j] - mean->data[0][j];
        }
    }
    
    Matrix* X_centered_t = matrix_transpose(X_centered);
    Matrix* cov = matrix_multiply(X_centered_t, X_centered);
    
    for (int i = 0; i < n_features; i++) {
        for (int j = 0; j < n_features; j++) {
            cov->data[i][j] /= (n_samples - 1);
        }
    }
    
    matrix_destroy(mean);
    matrix_destroy(X_centered);
    matrix_destroy(X_centered_t);
    
    return cov;
}

PCA* pca_fit(Matrix* X, int n_components) {
    if (!X || n_components <= 0 || n_components > X->cols) return NULL;
    
    PCA* pca = (PCA*)malloc(sizeof(PCA));
    if (!pca) return NULL;
    
    int n_samples = X->rows;
    int n_features = X->cols;
    
    pca->mean = matrix_create(1, n_features);
    for (int j = 0; j < n_features; j++) {
        double sum = 0.0;
        for (int i = 0; i < n_samples; i++) {
            sum += X->data[i][j];
        }
        pca->mean->data[0][j] = sum / n_samples;
    }
    
    Matrix* X_centered = matrix_create(n_samples, n_features);
    for (int i = 0; i < n_samples; i++) {
        for (int j = 0; j < n_features; j++) {
            X_centered->data[i][j] = X->data[i][j] - pca->mean->data[0][j];
        }
    }
    
    SVD* svd = svd_compute(X_centered);
    
    pca->n_components = n_components;
    pca->components = matrix_create(n_components, n_features);
    pca->eigenvalues = (double*)malloc(n_components * sizeof(double));
    pca->explained_variance_ratio = (double*)malloc(n_components * sizeof(double));
    
    double total_variance = 0.0;
    for (int i = 0; i < svd->S->rows && i < svd->S->cols; i++) {
        double s = svd->S->data[i][i];
        total_variance += s * s / (n_samples - 1);
    }
    
    for (int i = 0; i < n_components; i++) {
        double s = svd->S->data[i][i];
        pca->eigenvalues[i] = s * s / (n_samples - 1);
        pca->explained_variance_ratio[i] = pca->eigenvalues[i] / total_variance;
        
        for (int j = 0; j < n_features; j++) {
            pca->components->data[i][j] = svd->V->data[j][i];
        }
    }
    
    matrix_destroy(X_centered);
    svd_destroy(svd);
    
    return pca;
}

Matrix* pca_transform(PCA* pca, Matrix* X) {
    if (!pca || !X || X->cols != pca->mean->cols) return NULL;
    
    int n_samples = X->rows;
    int n_features = X->cols;
    
    Matrix* X_centered = matrix_create(n_samples, n_features);
    for (int i = 0; i < n_samples; i++) {
        for (int j = 0; j < n_features; j++) {
            X_centered->data[i][j] = X->data[i][j] - pca->mean->data[0][j];
        }
    }
    
    Matrix* components_t = matrix_transpose(pca->components);
    Matrix* X_transformed = matrix_multiply(X_centered, components_t);
    
    matrix_destroy(X_centered);
    matrix_destroy(components_t);
    
    return X_transformed;
}

Matrix* pca_inverse_transform(PCA* pca, Matrix* X_transformed) {
    if (!pca || !X_transformed || X_transformed->cols != pca->n_components) return NULL;
    
    int n_samples = X_transformed->rows;
    int n_features = pca->components->cols;
    
    Matrix* X_reconstructed = matrix_multiply(X_transformed, pca->components);
    
    for (int i = 0; i < n_samples; i++) {
        for (int j = 0; j < n_features; j++) {
            X_reconstructed->data[i][j] += pca->mean->data[0][j];
        }
    }
    
    return X_reconstructed;
}

void pca_destroy(PCA* pca) {
    if (!pca) return;
    
    matrix_destroy(pca->components);
    matrix_destroy(pca->mean);
    free(pca->eigenvalues);
    free(pca->explained_variance_ratio);
    free(pca);
}