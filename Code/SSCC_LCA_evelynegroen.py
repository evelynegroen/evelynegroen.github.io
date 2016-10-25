# Procedure:     Global sensitivity analysis for matrix-based LCA
# Method:        Squared Spearman (rank) correlation coefficient (SSCC)
#                MCS: Monte Carlo simulation (normal random)
# Author:        Evelyne Groen {evelyne [dot] groen [at] gmail [dot] com}
# Last update:   25/10/2016


import numpy as np

# Deterministic (_det) result:
A_det = np.matrix('10 0; -2 100')       #A-matrix
B_det = np.matrix('1 10')               #B-matrix

f = np.matrix('1000; 0')                #Functional unit vector

g_LCA = B_det * A_det.I * f

print("Deterministic value:", g_LCA[0,0])

# Step 1: Monte Carlo simulation using normal distribution functions for all input parameters
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

#Step 2: Calculate the Spearman rank correlation coefficients
#Reshape the data

g_list = np.reshape([g[0,0] for g in gs], (N,1))
As_list = np.reshape(As, (N,4))
Bs_list = np.reshape(Bs, (N,2))

Ps_list = np.concatenate((np.ones((N,1)), As_list[:,:1], As_list[:,2:], Bs_list), axis=1)

from scipy.stats import rankdata
import statistics as stats

#Transform to rankdata
Ps_rank = [rankdata(Ps_list[:,k]) for k in range(1,6)]
g_rank = [rankdata(g_list[:,0])]

covar = [np.cov(Ps_rank[k],g_rank[0]) for k in range(0,5)]

std_rank = np.std(Ps_rank, axis = 1)
std_rank_g = np.std(g_rank)

SSCC = [(covar[k][0,1] / (std_rank[k] * std_rank_g))**2 for k in range(0,5)]
print("squared Spearman correlation coefficients:", SSCC)

#Visualize: make a bar plot
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt

SSCC_procent = [SSCC[k] * 100 for k in range(0,5)]

x_label=[ 'A(1,1)', 'A(2,1)', 'A(2,2)', 'B(1,1)', 'B(1,2)']
x_pos = range(5)
plt.bar(x_pos, SSCC_procent, align='center')
plt.xticks(x_pos, x_label)
plt.title('Global sensitivity analysis: squared Spearman correlation coefficients')
plt.ylabel('SSCC (%)')
plt.xlabel('Parameter')
plt.show('Figure 1')
