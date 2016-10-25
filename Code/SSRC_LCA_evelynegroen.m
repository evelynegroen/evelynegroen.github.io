%Procedure:     Global sensitivity analysis (GSA) for matrix-based LCA
%Method:        GSA: Squared standardized regression coefficients  % 
%               Uncertainty propagation: Monte Carlo simulation
%Author:        Evelyne Groen {evelyne [dot] groen [at] gmail [dot] com}
%Last update:   20/10/2016
%Toolbox:       statistics_toolbox


A_det=[10 0; -2 100];       %A-matrix
B_det=[1 10];               %B-matrix
f=[1000; 0];                %Scaling vector f
g_LCA=B_det*(A_det\f);      %Deterministic result  

CV=0.05;                    %CV: coefficient of variation

%Step 1: Uncertianty propagation: Monte Carlo simulation
N=1000;                     %Sample size

A=zeros(2,2,N);             %Pre-allocate zeros for speed
B=zeros(1,2,N);
g=zeros(N,1);


for i=1:N; 
        A(:,:,i)=normrnd(A_det, abs(CV*A_det)); %avoid negative variation for A
        B(:,:,i)=normrnd(B_det, abs(CV*B_det));
        g(i)=B(:,:,i)*(A(:,:,i)\f);
end

%Step 2: Global sensitivity analysis: squared standardized regression coefficients
A_res=zeros(numel(A_det),N);  %Pre-allocate zeros for speed
B_res=zeros(numel(B_det),N);

for i=1:N
    A_res(:,i)=reshape(A(:,:,i),1,numel(A_det));
    B_res(:,i)=reshape(B(:,:,i),1,numel(B_det));
end

A_res(~any(A_res,2),:) =[];     % Delete samples for parameters equal to zero in A before calcultating regression coefficients
                                % NB: If B also contains parameters equal
                                % to zero, remove them as well, otherwise
                                % it is not possible to calculate the
                                % inverse.
P=[ones(N,1) A_res' B_res'];    % Combine samples for each input parameter in one matrix P; add ones to perform linear regression

%regression coefficients (RC):

%RC = inv(P'*P)*P'*g;          % This is what "regress" is doing in the
                               % next line
RC = regress(g, P);            % Regression coefficients for each non-zero parameters, first element equals the intersect

SSRC_A=zeros(nnz(A_det),1);
SSRC_B=zeros(nnz(B_det),1);

for j=1:nnz(A_det);
    SSRC_A(j)=((RC(j+1))^2)*var(A_res(j,:))/var(g);
end

for j=1:numel(B_det);
    SSRC_B(j)=((RC(nnz(A_det)+j+1))^2)*var(B_res(j,:))/var(g);
end

elem_AB = [find(A_det); find(B_det')];
SSRC=[SSRC_A; SSRC_B];

SRC_squared=[SSRC elem_AB]

