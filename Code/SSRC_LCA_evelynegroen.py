# Procedure:     Uncertainty propagation for matrix-based LCA
# Method:        MCS: Monte Carlo simulation (normal random)
# Author:        Evelyne Groen {evelyne [dot] groen [at] gmail [dot] com}
# Last update:   12/10/2016

import numpy as np

A_det = np.matrix('10 0; -2 100')       #A-matrix
B_det = np.matrix('1 10')               #B-matrix

f = np.matrix('1000; 0')                #Functional unit vector

g_LCA = B_det * A_det.I * f

print(g_LCA[0,0])                       #Deterministic result

# Monte Carlo simulation using normal distribution functions for all input parameters
# The mean values are equal to the initial values of A and B.
# The standard deviation equals 5% of the mean of A and B.

N = 10                                #Sample size
CV = 0.05                               #Coefficient of variation (CV = sigma/abs(mu)),

import random

As = [random.normalvariate(A_det, CV*A_det) for i in range(N)]
Bs = [random.normalvariate(B_det, CV*B_det) for i in range(N)]
f  = np.matrix('1000; 0')

gs = [B * A.I * f for A, B in zip(As, Bs)]

g_list =[g[0,0] for g in gs]

import scipy as sc

As_list = np.reshape(As, (N,4))
Bs_list = np.reshape(Bs, (N,2))

Ps_list = np.concatenate((np.ones((N,1)), As_list, Bs_list), axis=1)
#print(type(Ps_list), Ps_list)

Ps_list_new = np.concatenate((Ps_list[:, 0:2], Ps_list[:, 3:7]), axis=1) ##tja....
print(Ps_list_new)
#RC = ( ( (Ps_list.T*Ps_list) ).I )*(P_list.T)*g_list  #RC = inv(P'*P)*P'*g;
RC = np.inverse( (np.dot(Ps_list.T, Ps_list)) )


print(RC)
#RC = sc.reg(g_list, P_list)
