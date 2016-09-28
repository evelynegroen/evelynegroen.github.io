%Procedure:     Global sensitivity analysis (GSA) for matrix-based LCA
%Method:        GSA: Squared Spearman correlation coefficients  % 
%               Uncertainty propagation: Monte Carlo simulation
%Author:        Evelyne Groen {evelyne [dot] groen [at] gmail [dot] com}
%Last update:   26/09/2016 


A_det=[10 0; -2 100];       %A-matrix
B_det=[1 10];               %B-matrix
f=[1000; 0];                %Scaling vector f
g_LCA=B_det*(A_det\f);      %Deterministic result  

CV=0.05;                    %CV: coefficient of variation

%Step 1: Uncertianty propagation: Monte Carlo simulation

N=1000;                      %Sample size

A=zeros(2,2,N);             %Pre-allocate zeros for speed
B=zeros(1,2,N);
g=zeros(N,1);


for i=1:N; 
        A(:,:,i)=normrnd(A_det, abs(CV*A_det)); %avoid negative variation for A
        B(:,:,i)=normrnd(B_det, abs(CV*B_det));
        g(i)=B(:,:,i)*(A(:,:,i)\f);
end

%Step 2: Global sensitivity analysis: squared Spearman correlation coefficients

A_res=zeros(numel(A_det),N);    %Pre-allocate zeros for speed
B_res=zeros(numel(B_det),N);

for i=1:N
    A_res(:,i)=reshape(A(:,:,i),1,numel(A_det));
    B_res(:,i)=reshape(B(:,:,i),1,numel(B_det));
end

P=[A_res; B_res];               %Combine samples for each input parameter in one matrix P

rho_Spearman=corr(P',g,'type', 'Spearman'); %NB: rank correlation %
rho_Spearman(isnan(rho_Spearman))=[];

%squared Spearman correlation coefficients:
elem_AB = [find(A_det); find(B_det')];
SCC_squared=[rho_Spearman.^2 elem_AB]

