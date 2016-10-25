%Procedure:     Global sensitivity analysis (GSA) for matrix-based LCA with
%               correlated input parameters
%Method:        Sampling GSA
%               Latin hypercube sampling, covariance structure
%Author:        Evelyne Groen {evelyne [dot] groen [at] gmail [dot] com}
%Last update:   20/10/2016 
%Toolbox:       statistics_toolbox

A=[10 0; -2 100];       %A-matrix
B=[1 10];               %B-matrix
f=[1000; 0];            %Functional unit vector f

g=B*(A\f);

% Electricity production example Reinout %

A=[10 0; -2 100]; % 1 3; 2 4%
B=[1 10];          % 5 6%
f=[1000; 0];
g=B*(A\f);

stdA=[0.12 0; 0.23 12]; % standaardafwijkingen: Deze kun je aanpassen!%
varA=stdA.^2;
stdB=[0.01 1.2]; % standaardafwijkingen: Deze kun je aanpassen!%
varB=stdB.^2;
rho1=0.8; % Correlatie tussen A(1) en B(1): deze kun je ook aanpassen!
rho2=0; % Correlatie tussen A(4) en B(2): deze kun je ook aanpassen!
% rho1=1;
% rho2=1;

covar1=rho1*stdA(1)*stdB(1); %covariance between A(1) and B(1)
covar2=rho2*stdA(4)*stdB(2); %covariance between A(4) and B(2)


%Step 1: uncertainty propagation using a covariance matrix %%

N=5000; %sample size%
A_new=[10 -2 100]; %slordig...%
mu=[A_new(:);B(:)];
kT=nnz(stdA)+nnz(stdB);
corr=eye(kT);
corr(1,4)=rho1; corr(4,1)=rho1; % One less due to zero in A%
corr(3,5)=rho2; corr(5,3)=rho2;
stdA(stdA==0)=[];
std=[stdA(:); stdB(:)];
sigma=corr.*(std*std'); %Variance-covariance matrix
[T,num] = cholcov(sigma); % to test if lhsnorm can be used %

[X,P]=lhsnorm(mu,sigma,N); %Correlated sample design

zsA=find(A);
zsB=find(B);

A3=A; %Needed to put the correlated values back into the matrices
kA=nnz(A);
B3=B;
kB=nnz(B);
g_LHS=zeros(N,1);

for i=1:N; 
        A3(zsA)=X(i,1:kA);
        B3(zsB)=X(i,kA+1:kA+kB);
        g_LHS(i)=B3*(A3\f);
end

% TOTAL EXPLAINED VARIANCE BY EACH INPUT PARAMETER %

% Regression coefficients: regression y on only x(j), j= 1, 2, 3.. etc.
% j=parameter

for j=1:kT %number of parameters
    R(:,:,j)=[ones(N,1) X(:,j)];
    theta_hat(:,j)=regress(g_LHS,R(:,:,j)); %equation (12), Xu and Gertner, 2008
end

%Explained y by x(j)%
for j=1:kT
    y_hat(:,j)=theta_hat(1,j)+theta_hat(2,j)*X(:,j); %equation below equation (11)
end

%Total partial variance x(j)%
for j=1:kT
    V_hat(j)=(1/(N-1))*sum((y_hat(:,j)-mean(g_LHS)).^2); %equation (11)
end

% EXPLAINED VARIANCE BY UNCORRELATED VARIATION %

for j=1:kT 
    C=X;
    C(:,j)=[]; %remove column X(:,j)
    R_eta(:,:,j)=[ones(N,1) C];
    eta_hat(:,j)=regress(X(:,j),R_eta(:,:,j)); %equation (15) ??? fout?
end
   
%Estimated residual from regression over  all X(j) but j%
for j=1:kT
    C=X;
    C(:,j)=[]; %remove column X(:,j)
    z_hat(:,j)=X(:,j)-(eta_hat(1,j)+(C*eta_hat(2:end,j))); %equation (14)
end

%new regress uncorr%
for j=1:kT
    A_r(:,:,j)=[ones(N,1) z_hat(:,j)];
    r(:,j)=regress(g_LHS, A_r(:,:,j)); %equation (17)
end

%Explained variance of uncorrelated parameters%
for j=1:kT
    M(:,:,j)=[ones(N,1) z_hat(:,j)];
    y_hat_uncorr(:,j)=M(:,:,j)*r(:,j);
end

%Partial  variance of uncorrelated parameters%
for j=1:kT
    V_hat_uncorr(j)=(1/(N-1))*sum((y_hat_uncorr(:,j)-mean(g_LHS)).^2);
end

%Partial variance correlated parameters %
for j=1:kT
    V_hat_corr(j)=V_hat(j)-V_hat_uncorr(j);
end

%Sensitivity indices%
for j=1:kT
    S_LHS(j)=V_hat(j)/var(g_LHS);
    S_uncorr_LHS(j)=V_hat_uncorr(j)/var(g_LHS);
    S_corr_LHS(j)=V_hat_corr(j)/var(g_LHS);
end

% [L_S,IL_S]=sort(S);
% [L_Suncorr,IL_Suncorr]=sort(S_uncorr_LHS);
% [L_corr,ILcorr]=sort(abs(S_corr_LHS));

z2=[S_LHS;S_uncorr_LHS;S_corr_LHS]';
z4=[0 0 0]
z3=[z2(1:2,:);z4;z2(3:5,:)]

figure(2)
bar(z3)
title('Sensitivity indices for numerical approach using LHS')
legend('S_{tot}','S_{uncorr}','S_{corr}')
xlabel('parameter')
ylabel('percent/100')


