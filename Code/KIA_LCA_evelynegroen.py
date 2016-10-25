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

# Step Extra: Calculate the multipliers (You don't need this for KIA)

Gamma_A = np.multiply((A_det/g_LCA), dgdA)      # Multipliers of the A-matrix (for free)
print("The multipliers of the A-matrix are:")
print(Gamma_A)

Gamma_B = np.multiply((B_det/g_LCA), dgdB)      # Multipliers of the B-matrix (for free)
print("The multipliers of the B-matrix are")
print(Gamma_B)

# Step 2: Determine output variance (Just like AUP)

CV = 0.05                               # Coefficient of variation set to 5% (CV = sigma/mu)
                                        # You can easily change the CV to vary between 0 - 30%
                                        # For values CV > 30% you might get outliers

var_A = np.power(abs(CV*A_det), 2)      # Variance of the A-matrix (var =sigma^2),
var_B = np.power(abs(CV*B_det), 2)      # Variance of the B-matrix

P = np.concatenate((np.reshape(dgdA, 4), dgdB), axis=1)         # P contains partial derivatives of both A and B
var_P = np.concatenate((np.reshape(var_A, 4), var_B), axis=1)   # var_P contains all variances of each parameter in A and B

var_g = np.sum(np.multiply(np.power(P, 2), var_P))              # Total output variance (first order Taylor)

print("The output variance equals:", var_g)

# Step 3: Key issue analysis (variance decomposition)"
KIA = np.multiply(np.power(P, 2), var_P)/var_g
KIA_procent = [KIA[0][0,k]*100 for k in range(0,6)]
print("The contribution to the output variance of each intput parameter equals (in %):")
print(KIA_procent)

#Visualize: make a bar plot
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt

x_label=[ 'A(1,1)', 'A(1,2)', 'A(2,1)', 'A(2,2)', 'B(1,1)', 'B(1,2)']
x_pos = range(6)
plt.bar(x_pos, KIA_procent, align='center')
plt.xticks(x_pos, x_label)
plt.title('Global sensitivity analysis: squared standardized regression coefficients')
plt.ylabel('KIA (%)')
plt.xlabel('Parameter')
plt.show()
