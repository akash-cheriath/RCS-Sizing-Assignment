// Defines the two-phase simplex solver (linear programming)

// As the thrust >= 0 with a limit on maximum value, linear programming is required to select the thruster and thr duration of its burn
// to produce a required force and torque (simple matrix inversion is not sufficient).
// This solver helps in assesssing the control authority (largest pure force/torque per axis, with any 1 or 2 thrusters failed) helping
// in layout trade.

// Structure : minimize c^T x, subject to A x = b, x >= 0

#pragma once
#include <vector>
#include <cmath>

struct SimplexResult {
    bool ok = false;        // True when an optimum is found
    std::vector<double> x;  // holds optimum solution
    double obj = 0;         // holds optimum objective value
};

inline SimplexResult solveSimplex(const std::vector<std::vector<double>>& A,
                                  const std::vector<double>& b,
                                  const std::vector<double>& c) {

    const double eps = 1e-9;
    const int m = static_cast<int>(A.size());
    const int n = static_cast<int>(c.size());
    const int N = m + n;    // Total number of columns in the simplex table (Number of slack variable (= inequality equations) + variables)

    // In the simplex table rows 0 to m-1 are constraints and row m is the reduced-cost row
    // Column N holds the right-hand side value of the equations
    std::vector<std::vector<double>> T(m + 1, std::vector<double>(N + 1, 0.0));
    std::vector<int> basis(m);

    for (int i = 0; i < m; ++i) {
        double s = (b[i] < 0) ? -1.0 : 1.0 ; // To ensure all the RHS non-negative
        for (int j = 0; j < n; ++j)
            T[i][j] = s * A[i][j];
        T[i][n + i] = 1.0;
        T[i][N] = s * b[i];
        basis[i] = n + i;
    }

    auto pivot = [&](int r, int col) {
        double p = T[r][col];
        for (double& v : T[r])
            v /= p;
        for (int i = 0; i <= m; ++i) {
            if (i == r || std::fabs(T[i][col]) < 1e-15)
                continue;
            double f = T[i][col];
            for (int j = 0; j <= N; ++j)
                T[i][j] -= f * T[r][j];
        }
        basis[r] = col;
    };

    // Simplex iterations ; only columns < ncol may enter the basis
    auto iterate = [&](int ncol) -> int {
        for (int it = 0; it < 5000; ++it) {
            int col = -1;
            for (int j = 0; j < ncol; ++j) {    // Bland: first improving column
                if (T[m][j] < -eps) {
                    col = j;
                    break;
                }
            }
            if (col < 0)
                return 0;                       // Optimal
            int r = -1;
            double best = 1e300;
            for (int i = 0; i < m; ++i) {
                if (T[i][col] > eps) {
                    double ratio = T[i][N] / T[i][col];
                    bool tie = (r >= 0) && std::fabs(ratio - best) <= 1e-12 && basis[i] < basis[r];
                    if (ratio < best - 1e-12 || tie) {
                        best = ratio;
                        r = i;
                    }
                }
            }
            if (r < 0)
                return 1;                       // unbounded
            pivot(r, col);
        }
        return 2;                               // iteration limit
    };

    SimplexResult res;
    // PHASE 1 : Minimising the sum of artificaial variables
    for (int j = 0; j <= N; ++j) {
        if (j >= n && j < N) {
            T[m][j] = 0.0;
            continue;
        }
        double s = 0;
        for (int i =0; i < m; ++i)
            s += T[i][j];
        T[m][j] = -s;
    }

    iterate(N);

    // Feasibility test RELATIVE to the size of the right-hand side: with torque demands of ~1e3 Nm, round off 
    // leaves residuals far above any absolute 1e-7
    double bScale = 1.0;
    for (double v : b)
        bScale += std::fabs(v);
    
    if (-T[m][N] > 1e-9 * bScale)
        return res;                             // Infeasible
    
    // Driving remaining artificial variables out of the basis
    for (int i = 0; i < m; ++i) {
        if (basis[i] < n)
            continue;
        for (int j = 0; j < n; ++j) {
            if (std::fabs(T[i][j]) > eps) {
                pivot(i, j);
                break;
            }
        }
        // If no column found, the row is redundant; the artificial slack variable stays at 0
    }

    // PHASE 2 : Original objective, artificial columns should not reenter
    for (int j = 0; j <= N; ++j)
        T[m][j] = (j < n) ? c[j] : 0.0;
    for (int i = 0; i < m; ++i) {
        int k = basis[i];
        if (k < n && std::fabs(c[k]) > 0) {
            for (int j = 0; j <= N; ++j) {
                T[m][j] -= c[k] * T[i][j];
            }
        }
    }

    int st = iterate(n);
    if (st != 0)
        return res;                             // Unbounded or iteration limit
    
    res.x.assign(n, 0.0);
    for (int i = 0; i < m; ++i) {
        if (basis[i] < n)
            res.x[basis[i]] = T[i][N];
    }
    
    res.obj = 0;
    for (int j = 0; j < n; ++j)
        res.obj += c[j] * res.x[j];
    
    res.ok = true;

    return res;
}

