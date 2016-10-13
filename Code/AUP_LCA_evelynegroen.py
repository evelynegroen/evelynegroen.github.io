# Procedure:     Uncertainty propagation for matrix-based LCA
# Method:        AUP: Analytical Uncertainty Propagation (Taylor approximation, 1st order)
# Author:        Evelyne Groen {evelyne [dot] groen [at] gmail [dot] com}
# Last update:   12/10/2016

import numpy as np

A_det = np.matrix('10 0; -2 100')       # A-matrix
B_det = np.matrix('1 10')               # B-matrix

f = np.matrix('1000; 0')                # Functional unit vector

g_LCA = B_det * A_det.I * f

print("The deterministic result equals:", g_LCA[0,0])

# Step 1: Calculate partial derivatives
# N.B.: this is a vectorized implementation of the MatLab code that was
# originally written by Reinout Heijungs & Sangwong Suh

s = A_det.I * f                         # scaling vector s: inverse(A)*f
Lambda = B_det * A_det.I                # B*inverse(A)

dgdA = -(s * Lambda).T                  # Partial derivatives A-matrix
dgdB = s.T                              # Partial derivatives B-matrix

# Step 2: Determine output variance

CV = 0.05                               # Coefficient of variation set to 5% (CV = sigma/mu)
                                        # You can easily change the CV to vary between 0 - 30%
                                        # For values CV > 30% you might get outliers

var_A = np.power(abs(CV*A_det), 2)      # Variance of the A-matrix (var =sigma^2),
var_B = np.power(abs(CV*B_det), 2)      # Variance of the B-matrix

P = np.concatenate((np.reshape(dgdA, 4), dgdB), axis=1)         # P contains partial derivatives of both A and B
var_P = np.concatenate((np.reshape(var_A, 4), var_B), axis=1)   # var_P contains all variances of each parameter in A and B

var_g = np.sum(np.multiply(np.power(P, 2), var_P))              # Total output variance (first order Taylor)

print("The output variance equals:", var_g)
