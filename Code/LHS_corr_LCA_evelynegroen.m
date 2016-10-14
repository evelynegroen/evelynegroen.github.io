%Procedure:     Uncertainty propagation for matrix-based LCA with
%               correlated input parameters
%Method:        Latin hypercube sampling with covariance structure
%Author:        Evelyne Groen {evelyne [dot] groen [at] gmail [dot] com}
%Last update:   14/10/2016 
%Toolbox:       statistics_toolbox

A=[10 0; -2 100];       %A-matrix
B=[1 10];               %B-matrix
f=[1000; 0];            %Functional unit vector f

g=B*(A\f);

stdA=[0.12 0; 0.23 12]; % Standard distributions A-matrix%
stdB=[0.01 1.2];        

rho1=0.8;               % Correlation between A(1) and B(1)
rho2=0.9;               % Correlation between A(4) and B(2)

covar1=rho1*stdA(1)*stdB(1); % covariance between A(1) and B(1)
covar2=rho2*stdA(4)*stdB(2); % covariance between A(4) and B(2)

%LHS Correlated sample %%

N=5000;                         % Sample size
A_mu=[10 -2 100];               % Remove zero means%
mu=[A_mu(:);B(:)];

kT=nnz(stdA)+nnz(stdB);         % Total number of parameters with non-zero means
corr=eye(kT);
corr(1,4)=rho1; corr(4,1)=rho1; % NB "corr" contains only five nonzero parameters%
corr(3,5)=rho2; corr(5,3)=rho2;

stdA(stdA==0)=[];
stddev=[stdA(:); stdB(:)];

sigma=corr.*(stddev*stddev');   % Variance-covariance matrix
[T,num] = cholcov(sigma);       % to test if lhsnorm can be used 

[X,P]=lhsnorm(mu,sigma,N);      % Correlated sample design

zsA=find(A);
zsB=find(B);

A3=A;           %Needed to put the correlated values back into the matrices
kA=nnz(A);
B3=B;
kB=nnz(B);
g_LHS=zeros(N,1);

for i=1:N; 
        A3(zsA)=X(i,1:kA);
        B3(zsB)=X(i,kA+1:kA+kB);
        g_LHS(i)=B3*(A3\f);
end

meang=mean(g_LHS);      %mean
stdg=std(g_LHS);        %standard deviation