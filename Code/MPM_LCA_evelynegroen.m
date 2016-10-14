%Procedure:     Local sensitivity analysis for matrix-based LCA
%Method:        Multiplier method (MPM): first order Taylor (analytic)
%Author:        Evelyne Groen {evelyne [dot] groen [at] gmail [dot] com}
%Last update:   14/10/2016 
%Toolbox:       none

%NB:            This code is based on Heijungs and Suh (2002); http://dx.doi.org/10.1007/978-94-015-9900-9


A=[10 0; -2 100];   %A-matrix
B=[1 10];           %B-matrix
f=[1000; 0];        %Functional unit vector f
g_LCA=B*(A\f);      %Deterministic answer

s=A\f;              %Scaling vector s

%NB: this is a vectorized implementation of the original work of Reinout Heijungs & Sangwong Suh

s=A_det\f;                      %inv(A_det)*f
Lambda=B_det/A_det;             %B_det*inv(A)

dgdA=-Lambda'*s';               %Partial derivatives A-matrix
Gamma_A=(A_det./g_LCA).*dgdA    %Multipliers of the A-matrix

dgdB=s';                        %Partial derivatives B-matrix
Gamma_B=(B_det./g_LCA).*dgdB    %Multipliers of the B-matrix
