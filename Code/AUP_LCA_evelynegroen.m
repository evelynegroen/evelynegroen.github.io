%Procedure:     Uncertainty propagation for matrix-based LCA
%Method:        Analytic uncertainty propagation (AUP): first order Taylor (analytic)
%Author:        Evelyne Groen {evelyne [dot] groen [at] gmail [dot] com}
%Last update:   14/10/2016 
%Toolbox:       none

%NB:            This code is based on Heijungs and Suh (2002); http://dx.doi.org/10.1007/978-94-015-9900-9

A_det=[10 0; -2 100];       %A-matrix
B_det=[1 10];               %B-matrix
f=[1000; 0];                %Functional unit vector f
g_LCA=B_det*(A_det\f);      %Deterministic result  

CV=0.05;                    %CV: coefficient of variation (%)

var_A=abs(CV*A_det).^2;
var_B=abs(CV*B_det).^2;

%Step 1: Calculate partial derivatives
%NB: this is a vectorized implementation of the original work of Reinout Heijungs & Sangwong Suh

s=A_det\f;                      %inv(A_det)*f
Lambda=B_det/A_det;             %B_det*inv(A)

dgdA=-Lambda'*s';               %Partial derivatives A-matrix
Gamma_A=(A_det./g_LCA).*dgdA;   %For free: the multipliers of the A-matrix

dgdB=s';                        %Partial derivatives B-matrix
Gamma_B=(B_det./g_LCA).*dgdB;   %For free t(o)(w)o: the multipliers of the B-matrix

%Step 2: Determine output variance
P=[dgdA(:); dgdB(:)];           %P contains partial derivatives of both A and B 
var_P=[var_A(:); var_B(:)];     %var_P contaings the variances of each parameter in A and B

var_g=sum((P.^2).*var_P)       %Total output variance (first order Taylor)
