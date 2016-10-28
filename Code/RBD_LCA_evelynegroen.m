%Procedure:     Global sensitivity analysis for matrix-based LCA
%Method:        Random balance design
%Author:        Evelyne Groen {evelyne [dot] groen [at] gmail [dot] com}
%Last update:   25/10/2016 
%Toolbox:       statistics_toolbox

% NB: I adapted the MatLab code from "Global sensitivity analysis: the primer" (DOI: 10.1002/9780470725184), 
% where the code was given to perform a GSA using random balance
% design for a simple test case, written by S. Tarantola.

A = [10 0; -2 100];   % A-matrix
B = [1 10];           % B-matrix
f = [1000; 0];        % Functional unit vector f

g_LCA = B*(A\f);      % Total CO_2 emissions (deterministic result)

%Step 1: Draw equi-distributed values from normal distribution function 

N = 1000;                       % Sample size%
k = numel(A)+numel(B);          % number of parameters%
CV = 0.05;                      % Coefficient of variance
mu = [A(:)' B(:)'];

s0 = -pi:2*pi/(N-1):pi;         % Divide the input space of the frequency 
                                % space in the amount of samples 

s = zeros(N,k);                 % Pre-allocate for speed
x = zeros(N,k);                  
x_N = zeros(N,k);

for i=1:k; 
    s(:,i)=s0(randperm(N));     % Performs a random permutation for each parameter 
end

for i=1:k;                      % Sample from uniform distrubtion %
    x(:,i) = 0.5+asin(sin(s(:,i)))/pi; 
end

for i=1:k                       % Convert to normal distribution %
    x_N(:,i)=norminv(x(:,i),mu(i),abs(CV*mu(i))); 
end

x_N(isnan(x_N)) = 0;            % Replace NaN by 0 for parameters equal to zero
x_A = x_N(:,1:numel(A));        % Return to one for matrix A and one for matrix B
x_B = x_N(:,numel(A)+1:k);

%Step 2: perform uncertainty propagation

k_A = 1:numel(A);               % Pre-allocate for speed
k_B = 1:numel(B);
A_temp = zeros(size(A));
B_temp = zeros(size(B));
gs = zeros(size(g_LCA));

for j=1:N;
        A_temp(k_A) = x_A(j,:); % The non-zero elements of A are replaced random draws      
        B_temp(k_B) = x_B(j,:); % The same is done for the B-matrix %              
        g_temp = B_temp*(A_temp\f);
        gs(:,:,j) = g_temp;
end

%Step 3: Global sensitivity analysis

index = zeros(N,k);             % Pre-allocate for speed
dummy = zeros(N,k);
g_sort = zeros(N,k);
spectrum = zeros(N,k);
V = zeros(1,k);
S = zeros(1,k);

for i=1:k;                      % Orders the elements of s in ascending order
    [dummy(:,i),index(:,i)] = sort(s(:,i)); 
end
    
g_row = 1;                      % Choose emisison type, equal 1 if B contains only one row

for i=1:k;
    g_sort(:,i) = gs(g_row,(index(:,i)));             
end

for i=1:k;                      % fft is the fast Fourier transform
    spectrum(:,i) = (abs(fft(g_sort(:,i)))).^2/N; 
end

M = 6;                          % Harmonics 

for i=1:k;                      % See equation 4.28 in Global Sensitivity Analysis: the primer
    V(i)=2*sum(spectrum(2:M+1,i)); 
    Vtot=sum(spectrum(2:N));
    S(i)=V(i)/Vtot;
end

