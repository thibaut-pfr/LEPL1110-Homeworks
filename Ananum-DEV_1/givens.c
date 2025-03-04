    // Apply Givens rotations to zero out the sub-diagonals
    for (int j = k; j < n - 1; j++) {
        for (int i = j; i > j - k + 1; i--) {
            double a = get(A,j,j);
            double b = get(A,i,j);
            double r = sqrt(a * a + b * b);
            double c = a / r;
            double s = -b / r;
            // Apply rotation to rows i and j
            double start = MAX(k,i);
            double stop = MIN(i + k, n);
            for (int l = start; l < stop; l++) {
                double temp_i = l < j ? c * get(A,i,l) : c * get(A,i,l) - s * get(A,j,l);
                set(A,i,l, temp_i);
            }

            start = MAX(k,j);
            stop = MIN(j + k, n);
            for (int l = start; l < stop; l++) {
                double temp_j = l > i + k ? c * get(A,j,l) : s * get(A,i,l) + c * get(A,j,l);
                set(A,j,l, temp_j);
            }

            // Apply rotation to columns i and j
            start = MAX(0,i - k);
            stop = MIN(i + k, n);
            for (int l = start; l < stop; l++) {
                double temp_i = l < j ? c * get(A,l,i) : c * get(A,l,i) - s * get(A,l,j);
                set(A,l,i, temp_i);
            }

            start = MAX(0,j - k);
            stop = MIN(j + k, n);
            for (int l = start; l < stop; l++) {
                double temp_j = l > i + k ? c * get(A,l,j) : s * get(A,l,i) + c * get(A,l,j);
                set(A,l,j, temp_j);
            }

        }
    }