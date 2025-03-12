import numpy as np
from scipy.linalg import hessenberg
# Example usage:
A = np.array([[4, 1, 2, 3],
              [3, 4, 1, 2],
              [2, 3, 4, 1],
              [1, 2, 3, 4]], dtype=float)

# Example usage for band-symmetric matrix:
A_band = np.array([[10, 10, 10, 10, 10, 10, 10, 0, 0, 0],
                   [10, 9, 9, 9, 9, 9, 9, 9, 0, 0],
                   [10, 9, 8, 8, 8, 8, 8, 8, 8, 0],
                   [10, 9, 8, 7, 7, 7, 7, 7, 7, 7],
                   [10, 9, 8, 7, 6, 6, 6, 6, 6, 6],
                   [10, 9, 8, 7, 6, 5, 5, 5, 5, 5],
                   [10, 9, 8, 7, 6, 5, 4, 4, 4, 4],
                   [0, 9, 8, 7, 6, 5, 4, 3, 3, 3],
                   [0, 0, 8, 7, 6, 5, 4, 3, 2, 2],
                   [0, 0, 0, 7, 6, 5, 4, 3, 2, 1]], dtype=float)

A_band_small = np.array([[8,1,4,0],
                         [1,4,3,1],
                         [4,3,1,2],
                         [0,1,2,9]], dtype=float)

MATRIX = A_band

eigenvalues = np.linalg.eigvals(MATRIX)
print("Eigenvalues of H_band:")
print(eigenvalues)


# Read the true eigenvalues from the file
with open('hessenberg_matrix.txt', 'r') as file:
    linesplitted = [line.strip().split() for line in file]
    hessenberg_computed = np.array([[float(val) for val in line] for line in linesplitted])
    true_eigenvalues = np.linalg.eigvals(hessenberg_computed)
    print("Hessenberg matrix computed (rounded to 2 digits):")
    print(np.round(hessenberg_computed, 2))
    print("True eigenvalues:")
    print(true_eigenvalues)

    # Check if the computed eigenvalues match the true ones
if np.allclose(np.sort(eigenvalues), np.sort(true_eigenvalues), atol=1e-6):
    print("The computed eigenvalues match the true eigenvalues !")
else:
    print("The computed eigenvalues do not match the true eigenvalues.")
