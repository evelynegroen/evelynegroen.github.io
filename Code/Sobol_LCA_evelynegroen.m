%Procedure:     Global sensitivity analysis for matrix-based LCA
%Method:        Sobol' indicues
%               For the uncertainty propagation: Monte Carlo simulation (normal random)
%Author:        Evelyne Groen {evelyne [dot] groen [at] gmail [dot] com}
%Last update:   26/10/2016 
%Toolbox:       statistics_toolbox
    
A=[10 0; -2 100];                   % A-matrix
B=[1 10];                           % B-matrix
f=[1000; 0];                        % Functional unit vector f

g_LCA = B*(A\f);                    % Determnistic result

k = numel(A)+numel(B);              % Number of parameters
N = 1000;                           % Sample zie
CV = 0.05;                          % Coefficient of variation

%Step 1: create 3 matrices containing random draws for each parameter in A
%and B

%Step 1a: create two matrices containing independent samples %

A_rand_1 = zeros(N,numel(A));      % Pre-allocate for speed
A_rand_2 = zeros(N,numel(A));
B_rand_1 = zeros(N,numel(B));
B_rand_2 = zeros(N,numel(B));

for i=1:N; 
    for j=1:numel(A)
        A_rand_1(i,j) = normrnd(A(j), abs(CV*A(j)) );
        A_rand_2(i,j) = normrnd(A(j), abs(CV*A(j)) );     
    end
    for j=1:numel(B)
        B_rand_1(i,j) = normrnd(B(j), abs(CV*B(j)) );
        B_rand_2(i,j) = normrnd(B(j), abs(CV*B(j)) );
    end
end


AB_rand_1=[A_rand_1 B_rand_1];      % Sample matrix 1
AB_rand_2=[A_rand_2 B_rand_2];      % Sample matrix 2

% Step 1b: Create a third matrix based on sample matrix 1 & 2

D_rand = zeros(N,k,k);              % Pre-allocate for speed

for j=1:k;
    D_rand(:,:,j) = AB_rand_1;      % Create k  matrices filled with sample matrix 1
end

for j=1:k;
    D_rand(:,j,j) = AB_rand_2(:,j); % Replace the jth column with sample matrix 2
end

%D1=zeros(N,numel(A),k);
%D2=zeros(N,numel(B),k);

D_A = D_rand(:,1:numel(A),:) ;       % Split up into matrices A and B for LCA%
D_B = D_rand(:,numel(A)+1:numel(A)+numel(B),:);

%Step 2: peform uncertainty propagation with the 3 sample matrices

k_A = 1:numel(A);                   % Number of parameters in A
k_B = 1:numel(B);                   % Number of parameters in B

%Step 2a:  Calculate model output for sample matrix 1%

A_temp_1 = zeros(2,2);              % Pre-allocate for speed
B_temp_1 = zeros(1,2);
A_store_1 = zeros(2,2,N);
B_store_1 = zeros(1,2,N);
g_rand_1 = zeros(1,N);

for i=1:N;                  
    A_temp_1(k_A) = A_rand_1(i,:);
    A_store_1(:,:,i) = A_temp_1;
    
    B_temp_1(k_B) = B_rand_1(i,:);
    B_store_1(:,:,i) = B_temp_1;
    
    g_rand_1(i) = B_store_1(:,:,i) * (A_store_1(:,:,i)\f);
end

%Step 2b:  Calculate model output for sample matrix 2%

A_temp_2 = zeros(2,2);              % Pre-allocate for speed
B_temp_2 = zeros(1,2);
A_store_2 = zeros(2,2,N);
B_store_2 = zeros(1,2,N);
g_rand_2=zeros(1,N);

for i=1:N;                 
    A_temp_2(k_A) = A_rand_2(i,:);
    A_store_2(:,:,i) = A_temp_2;
    
    B_temp_2(k_B) = B_rand_2(i,:); 
    B_store_2(:,:,i) = B_temp_2;
    
    g_rand_2(i)=B_store_2(:,:,i) * (A_store_2(:,:,i)\f);
end

%Step 2c:  Calculate model output for sample matrix 3%
A_temp_D = zeros(2,2);              % Pre-allocate for speed
B_temp_D = zeros(1,2);
A_store_D = zeros(2,2,N,k);
B_store_D = zeros(1,2,N,k);
g_rand_D = zeros(N,k);

for j=1:k
    for i=1:N                     
        A_temp_D(k_A) = D_A(i,:,j);
        A_store_D(:,:,i,j) = A_temp_D;
        
        B_temp_D(k_B) = D_B(i,:,j);
        B_store_D(:,:,i,j) = B_temp_D;
    
        g_rand_D(i,j) = B_store_D(:,:,i,j) * (A_store_D(:,:,i,j)\f);
    end
end

% Step 3: Calculate output variance

% Some statistics (mean and variance): 
mean_sobol_square = (1/N^2)*sum(g_rand_2).*sum(g_rand_1);
mean_sobol = sqrt(mean_sobol_square); 

var_g = (1/N)*sum(g_rand_1.*g_rand_1) - (1/N^2)*sum(g_rand_1).*sum(g_rand_1);                       %Total variance%

% Main sensitivity index (or first order sensitivity index) and total sensitivity index
S_main = zeros(1,k);            % Pre-allocate for speed
S_total =zeros(1,k);
var_p = zeros(1,k);             % "Expected reduction in variance that would be obtained if [parameter] k could be fixed" (Saltelli et al., 2010)
exp_p = zeros(1,k);             % "Expected variance that would be left if all [parameters] but [parameter] k could be fixed" (Saltelli et al., 2010)


for j=1:k;
    var_p(j)=var_g - (1/(2*N)) * sum((g_rand_2' - g_rand_D(:,j)).^2); % There are different versions for the Sobol' main effect, I found that this one worked best
    
    exp_p(j)=(1/(2*N))*sum((g_rand_1'-g_rand_D(:,j)).^2); 
    
    S_main(j) = var_p(j)/var_g;     % Main sensitivity index
    S_total(j) = exp_p(j)/var_g;    % Total sensitivity index
end



