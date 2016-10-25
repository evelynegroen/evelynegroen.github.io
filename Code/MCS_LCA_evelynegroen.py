# Procedure:     Uncertainty propagation for matrix-based LCA
# Method:        MCS: Monte Carlo simulation (normal random)
# Author:        Evelyne Groen {evelyne [dot] groen [at] gmail [dot] com}
# Last update:   25/10/2016

import numpy as np

#Deterministic (_det) result:
A_det = np.matrix('10 0; -2 100')       #A-matrix
B_det = np.matrix('1 10')               #B-matrix

f = np.matrix('1000; 0')                #Functional unit vector

g_LCA = B_det * A_det.I * f

print(g_LCA[0,0])                       #Deterministic result

# Monte Carlo simulation using normal distribution functions for all input parameters
# The mean values are equal to the initial values of A and B.
# The standard deviation equals 5% of the mean of A and B.

N = 1000                                #Sample size
CV = 0.05                               #Coefficient of variation (CV = sigma/abs(mu)),

import random


A1 = [random.gauss(A_det[0,0], CV*A_det[0,0]) for i in range(N)]
A3 = [random.gauss(A_det[1,0], CV*A_det[1,0]) for i in range(N)]
A4 = [random.gauss(A_det[1,1], CV*A_det[1,1]) for i in range(N)]

B1 = [random.gauss(B_det[0,0], CV*B_det[0,0]) for i in range(N)]
B2 = [random.gauss(B_det[0,1], CV*B_det[0,1]) for i in range(N)]


As = [np.matrix([[A1[i], 0],[A3[i], A4[i]]]) for i in range(N)]
Bs = [np.matrix([[B1[i], B2[i]]]) for i in range(N)]

f  = np.matrix('1000; 0')

gs = [B * A.I * f for A, B in zip(As, Bs)]

g_list =[g[0,0] for g in gs]

import statistics as stats
var_g = stats.variance(g_list)

print("output variance:", var_g)

import matplotlib
import matplotlib.pyplot as plt

plt.hist(g_list,20)

plt.title("Histogram")
plt.xlabel("kg CO2")
plt.ylabel("Frequency")

plt.show()
