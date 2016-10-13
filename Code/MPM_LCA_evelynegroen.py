# Procedure:     Local sensitivity analysis for matrix-based LCA
# Method:        MPM: Multiplier Method (Taylor approximation, 1st order)
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

s = A_det.I * f                                 # scaling vector s: inverse(A)*f
Lambda = B_det * A_det.I                        # B*inverse(A)

dgdA = -(s * Lambda).T                          # Partial derivatives A-matrix
dgdB = s.T                                      # Partial derivatives B-matrix

# Step 2: Calculate the multipliers

Gamma_A = np.multiply((A_det/g_LCA), dgdA)      # Multipliers of the A-matrix
print("The multipliers of the A-matrix are:")
print(Gamma_A)

Gamma_B = np.multiply((B_det/g_LCA), dgdB)      # Multipliers of the B-matrix
print("The multipliers of the B-matrix are")
print(Gamma_B)
