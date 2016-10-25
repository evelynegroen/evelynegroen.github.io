%Procedure:     Uncertainty propagation for matrix-based LCA
%Method:        Latin hypercube sampling (homemade)
%Author:        Evelyne Groen {evelyne [dot] groen [at] gmail [dot] com}
%Last update:   20/10/2016 
%Toolbox:       statistics_toolbox

%N.B: Alternative approach: use the build-in lhsnorm-function
%Open LHS_corr_LCA_evelynegroen.m and set the off-diagonal covariance to zero.

A=[10 0; -2 100];   %A-matrix
B=[1 10];           %B-matrix
f=[1000; 0];        %scaling vector f

g_LCA=B*(A\f);      %deterministic answer (kg CO_2 emissions)

n1=nnz(A);          %number of nonzero elements in A%
n2=nnz(B);          %number of nonzero elements in B%

[r,c,v]=find(A);    %r: row; c: column; v: value
[r2,c2,v2]=find(B);

CV=0.05;            %CV: coefficient of variation (set to 5%)

N=1000;             %Sample size %
D=zeros(n1,N);      %Pre-allocate for speed
E=zeros(n2,N);
Q=zeros(size(A,1),size(A,1));
P=zeros(size(B,1),2);
g=zeros(1,N);

for s=1:n1                  
        D(s,:)=randperm(N)-1;               % For each non-zero element in A a unique randperm is created equal to the amount of samples%
end
   
for t=1:n2 
        E(t,:)=randperm(N)-1;
end
    
for i=1:N
    x(:,i)=rand(n1,1)/N+D(:,i)/N;            % From each interval D/N a scaled probability is drawn %%
    Q1(:,:,i)=norminv(x(:,i),v, abs(CV*v));  % NB: it is not possible to take norminverse of sparse matrix %
    Q(A~=0)=Q1(:,:,i);                       % Create full matrix by replacing each non-zero element in Q (Q~=0) by the sampled value A %
    
    y(:,i)=rand(n2,1)/N+E(:,i)/N;
    P1(:,:,i)=norminv(y(:,i),v2',abs(CV*v2'));
    P(B~=0)=P1(:,:,i);
    
    g(i)=P*(Q\f);
end
    
meang=mean(g);
stdg=std(g);

figure(1)
hist(g)
hTitle=title('Latin hypercube sampling');
hXLabel=xlabel('kg CO_{2}e');
hYLabel=ylabel('Number of realizations');


