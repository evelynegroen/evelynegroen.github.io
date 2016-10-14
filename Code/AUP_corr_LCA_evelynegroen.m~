%Procedure:     Uncertainty propagation for matrix-based LCA with
%               correlated input parameters
%Method:        Analytic uncertainty propagation (AUP correlated): first order Taylor (analytic)
%Author:        Evelyne Groen {evelyne [dot] groen [at] gmail [dot] com}
%Last update:   14/10/2016 
%Toolbox:       none

A=[10 0; -2 100];       %A-matrix
B=[1 10];               %B-matrix
f=[1000; 0];            %Functional unit vector f

g=B*(A\f);

stdA=[0.12 0; 0.23 12]; % Standard distributions A-matrix%
varA=stdA.^2;           % Variance A-matrix
stdB=[0.01 1.2];        
varB=stdB.^2;
rho1=0.8;               % Correlation between A(1) and B(1)
rho2=0.9;               % Correlation between A(4) and B(2)

covar1=rho1*stdA(1)*stdB(1); % covariance between A(1) and B(1)
covar2=rho2*stdA(4)*stdB(2); % covariance between A(4) and B(2)

s=A\f;
Lambda=B*inv(A);

dgdA=-Lambda(:)*s';         % Partial derivatives A-matrix
dgdB=s';                    % Partial derivatives B-matrix


%Output variance ignoring correlations:
varg=sum(((dgdA(:)).^2).*varA(:))+sum(((dgdB(:)).^2).*varB(:))

%Output variance including correlations / Groen & Heijungs 2016: Equation (10)
varg_tot=sum(((dgdA(:)).^2).*varA(:))+sum(((dgdB(:)).^2).*varB(:))+2*dgdA(1)*dgdB(1)*covar1+2*dgdA(4)*dgdB(2)*covar2
