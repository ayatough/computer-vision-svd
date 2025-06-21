# SVD Applications in C

A comprehensive from-scratch implementation of Singular Value Decomposition (SVD) in C using the Jacobi eigenvalue algorithm, with multiple computer vision and linear algebra applications.

## Files

- `matrix.h` - Header file with all structures and function declarations
- `matrix.c` - Complete implementation (~918 lines)
- `test_svd.c` - Test program for basic SVD functionality
- `test_8point.c` - Test program for essential matrix estimation
- `test_homography.c` - Test program for homography algorithms
- `test_pca.c` - Test program for Principal Component Analysis
- `Makefile` - Build configuration with multiple targets

## Implemented Algorithms

### Core Matrix Operations
- Dynamic matrix allocation and deallocation
- Matrix transpose and multiplication
- Jacobi eigenvalue algorithm for symmetric matrices
- SVD decomposition: A = U × S × V^T
- Support for both square and rectangular matrices
- Matrix pseudoinverse computation
- 3×3 matrix determinant and inverse

### Computer Vision Applications

#### 1. Essential Matrix Estimation (8-Point Algorithm)
- Point correspondence data structures
- Coordinate normalization for numerical stability
- Constraint matrix construction from epipolar geometry
- SVD-based solution of homogeneous linear system
- Essential matrix constraint enforcement (rank 2, equal singular values)
- **Applications**: Stereo vision, camera calibration, 3D reconstruction

#### 2. Homography Estimation & Decomposition
- **DLT (Direct Linear Transform)** for 2D-2D homography estimation
- **Homography decomposition** into rotation, translation, and normal vectors
- Coordinate normalization and denormalization
- **Applications**: Image rectification, planar object tracking, augmented reality

### Mathematical Applications

#### 3. Principal Component Analysis (PCA)
- Covariance matrix computation
- SVD-based PCA fitting
- Data transformation and inverse transformation
- Explained variance ratio calculation
- **Applications**: Dimensionality reduction, feature extraction, data compression

## Building

```bash
make                  # Compile all programs
make run              # Run SVD tests
make run-8point       # Run 8-point algorithm tests
make run-homography   # Run homography tests
make run-pca          # Run PCA tests
make test-all         # Run all tests
make clean            # Clean build files
```

## Usage Examples

### SVD Computation
```c
#include "matrix.h"

// Create a matrix
Matrix* A = matrix_create(3, 2);
matrix_set(A, 0, 0, 1.0);
matrix_set(A, 0, 1, 2.0);
// ... set other values

// Compute SVD
SVD* svd = svd_compute(A);

// Access results
printf("Rank: %d\n", svd->rank);
matrix_print(svd->U);    // Left singular vectors
matrix_print(svd->S);    // Singular values (diagonal)
matrix_print(svd->V);    // Right singular vectors

// Cleanup
svd_destroy(svd);
matrix_destroy(A);
```

### Essential Matrix Estimation
```c
// Create point correspondences (need at least 8 points)
PointCorrespondences* corr = correspondences_create(8);

// Set point correspondences (x1, y1) <-> (x2, y2)
correspondences_set(corr, 0, 100, 200, 150, 180);
// ... set other correspondences

// Compute essential matrix using 8-point algorithm
Matrix* E = compute_essential_matrix(corr);

// Cleanup
correspondences_destroy(corr);
matrix_destroy(E);
```

### Homography Estimation
```c
// Create point correspondences (need at least 4 points)
PointCorrespondences* corr = correspondences_create(4);
// ... set correspondences

// Compute homography using DLT
Matrix* H = compute_homography_dlt(corr);

// Decompose homography (requires camera intrinsics)
Matrix* K1 = matrix_create(3, 3); // Camera intrinsics
// ... set K1
HomographyDecomposition* decomp = decompose_homography(H, K1, K1);

// Access rotation, translation, normal vectors
matrix_print(decomp->R1);
matrix_print(decomp->t1);
matrix_print(decomp->n1);

// Cleanup
homography_decomposition_destroy(decomp);
```

### Principal Component Analysis
```c
// Create data matrix (samples × features)
Matrix* X = matrix_create(100, 10);
// ... populate with data

// Fit PCA with 3 components
PCA* pca = pca_fit(X, 3);

// Transform data to reduced dimensions
Matrix* X_reduced = pca_transform(pca, X);

// Access explained variance
printf("Explained variance ratios: ");
for (int i = 0; i < pca->n_components; i++) {
    printf("%.3f ", pca->explained_variance_ratio[i]);
}

// Cleanup
pca_destroy(pca);
```

## Algorithm Details

### SVD Implementation
Uses eigendecomposition approach:
1. Compute A^T×A and A×A^T
2. Find eigenvalues/eigenvectors using Jacobi method
3. Singular values are square roots of eigenvalues
4. Construct U and V matrices from eigenvectors

### 8-Point Algorithm
For essential matrix estimation:
1. Normalize point coordinates for numerical stability
2. Build constraint matrix from epipolar constraint
3. Solve Ax = 0 using SVD (smallest singular vector)
4. Reshape to 3×3 matrix and enforce essential matrix constraints
5. Denormalize the result

### Homography DLT
For 2D-2D transformation estimation:
1. Normalize coordinates in both images
2. Build constraint matrix from point correspondences
3. Solve using SVD (smallest singular vector)
4. Denormalize resulting homography

### PCA via SVD
For dimensionality reduction:
1. Center the data (subtract mean)
2. Compute SVD of centered data matrix
3. Principal components are columns of V matrix
4. Eigenvalues are squares of singular values

## Test Results

All algorithms demonstrate:
- **Perfect reconstruction** for synthetic data
- **Numerical stability** through coordinate normalization
- **Correct constraint enforcement** (rank, orthogonality)
- **Accurate variance explained** in PCA
- **Zero projection error** for exact correspondences

## Limitations

- Uses double precision floating point
- Jacobi method convergence depends on tolerance (1e-12) and max iterations (1000)
- Memory allocation errors are handled but not extensively tested
- Performance not optimized for very large matrices
- No outlier rejection (RANSAC would be needed for real applications)
- Homography decomposition assumes planar scenes
- PCA assumes linear relationships in data

## Potential Extensions

- **Perspective-n-Point (PnP)** solver for camera pose estimation
- **Procrustes analysis** for point set alignment
- **Linear least squares** solver for overdetermined systems
- **RANSAC** for robust estimation with outliers
- **Bundle adjustment** for 3D reconstruction refinement