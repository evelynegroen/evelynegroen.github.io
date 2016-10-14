%Procedure:     Uncertainty propagation for matrix-based LCA
%Method:        Monte Carlo simulation (normal random)
%Author:        Evelyne Groen {evelyne [dot] groen [at] gmail [dot] com}
%Last update:   14/10/2016 
%Toolbox:       statistics_toolbox
 

A=[10 0; -2 100];   %A-matrix
B=[1 10];           %B-matrix
f=[1000; 0];        %Functional unit vector f

g_LCA=B*(A\f);      %Total CO_2 emissions (deterministic result)

N=1000;             %Sample size%

%Monte Carlo simulation assuming normal distributions (normrnd)
%The mean equals the initial values of A and B
%The standard deviation equals 5% of the mean of A and B

CV=0.05;            %Relative variation (CV=std/mu): 5%

for i=1:N; %repetitions %
    A_sam=normrnd(A, abs(CV*A));  %A_sam is the sampled matrix A
    B_sam=normrnd(B, abs(CV*B));  %B_sam is the sampled matrix B
    g(i)=B_sam*(A_sam\f);         %g(i) contains the result for each Nth iteration
end

meang=mean(g);      %mean
stdg=std(g);        %standard deviation

figure(1)
hist(g)
hTitle=title('Monte Carlo sampling');
hXLabel=xlabel('kg CO_{2}e');
hYLabel=ylabel('Number of realizations');

