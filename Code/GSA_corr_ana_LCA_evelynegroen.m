%Procedure:     Global sensitivity analysis (GSA) for matrix-based LCA with
%               correlated input parameters
%Method:        Analytic GSA
%               Analytic uncertainty propagation (AUP correlated): first order Taylor (analytic)
%Author:        Evelyne Groen {evelyne [dot] groen [at] gmail [dot] com}
%Last update:   20/10/2016 
%Toolbox:       none

A=[10 0; -2 100];       %A-matrix
B=[1 10];               %B-matrix
f=[1000; 0];            %Functional unit vector f

g=B*(A\f);

%Step 1: determine the input uncertainties

stdA=[0.12 0; 0.23 12]; % Standard deviations A-matrix%
varA=stdA.^2;           % Variance A-matrix
stdB=[0.01 1.2];        
varB=stdB.^2;
rho1=0.8;               % Correlation between A(1) and B(1)
rho2=0.9;               % Correlation between A(4) and B(2)

covar1=rho1*stdA(1)*stdB(1); % covariance between A(1) and B(1)
covar2=rho2*stdA(4)*stdB(2); % covariance between A(4) and B(2)

%Step 2: Determine the output variance
s=A\f;
Lambda=B*inv(A);

dgdA=-Lambda(:)*s';         % Partial derivatives A-matrix
dgdB=s';                    % Partial derivatives B-matrix


%Output variance ignoring correlations:
varg=sum(((dgdA(:)).^2).*varA(:))+sum(((dgdB(:)).^2).*varB(:));

%Output variance including correlations / Groen & Heijungs 2016: Equation (10)
varg_tot=sum(((dgdA(:)).^2).*varA(:))+sum(((dgdB(:)).^2).*varB(:))+2*dgdA(1)*dgdB(1)*covar1+2*dgdA(4)*dgdB(2)*covar2;

%Step 3: Calculate the sensitivity indices

%Step 3a: uncorrelated partial variance
dgdx=[dgdA(:); dgdB(:)];
rho=[rho1 0 0 rho2 rho1 rho2];
varx=[varA(:); varB(:)];

var_uncorr=(1-rho(:).^2).*((dgdx(:)).^2).*varx(:);

%Step 3b: correlated partial variance
var_corr1=2*dgdA(1)*dgdB(1)*covar1+rho1^2*((dgdA(1)).^2).*varA(1)+rho1^2*((dgdB(1)).^2).*varB(1);
var_corr2=2*dgdA(4)*dgdB(2)*covar2+rho2^2*((dgdA(4)).^2).*varA(4)+rho2^2*((dgdB(2)).^2).*varB(2);

var_corr=[var_corr1 0 0 var_corr2 var_corr1 var_corr2]';

%Step 3c: total partial variance

var_tot=var_corr+var_uncorr;

%Sensitivity indices%
S_tot=var_tot./varg_tot;
S_uncorr=var_uncorr./varg_tot;
S_corr=var_corr./varg_tot;

%Figure

z=[S_tot';S_uncorr';S_corr']';
bar(z)
title('Sensitivity indices for analytical approach')
legend('S_{tot}','S_{uncorr}','S_{corr}')
xlabel('parameter')
ylabel('percent/100')
