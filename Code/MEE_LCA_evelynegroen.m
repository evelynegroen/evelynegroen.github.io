%Method of elementary effects using matrix-based LCA % 
%Author: Evelyne Groen {evelyne [dot] groen [at] gmail [dot] com}

A=[10 0; -2 100];       %A-matrix %
B=[1 10];               %B-matrix %
f=[1000; 0];            %Functional unit vector f %
g_LCA=B*(A\f);          %Deterministic answer (kg CO_2) %

% Input for creating the trajectories %
k=nnz(A) +nnz(B);       % Amount of uncertain input parameters in A and B%
k1=k+1;                 % k+1, amout of steps through the grid in the hypercube %
z=10;                   % Amount of trajectories through the hypercube %
Ps=zeros(k1,k);         % Pre-allocate zeros to each variable for speed %
Ps_A=zeros(k1,nnz(A));
Ps_B=zeros(k1,nnz(B));

for t=1:z;
    % create a base for your trajectory, chosen at random %
    x1=randi([0,1],1,k);        % create an array of random integers between 0 and 1 of size k%
    x1(x1==1)=1/3;              % replace ones by 1/3, if you have a grid of 4 or 5, change this line of code and the line above %
    xi=diag(x1)  ;              % random sample of the base of the trajectory %
    
    % Create a random trajectory through the grid %
    D1=randi([0,1],[1,k]);      % create an array of random integers of size k between 0 and 1 %
    D1(~D1)=-1;                 % replace the zeros by -1%
    Ds=diag(D1) ;               % create a diagonal matrix which causes the trajectory to decrease or increase %

    PI=eye(k);                  % Identity matrix%
    Pr=PI(:,randperm(k));       % (random) gives order which of the input factors is changed with delta %

    B_grid=tril(ones(k+1,k),-1);    % size of the sampled grid (k+1)xk, makes sure that each parameter is analyzed %   
    J=ones(k1,k);                   % (k+1) x k size of the sample trajectory%

    delta=2/3;                      % size of the steps in the grid %
 
    Ps(:,:,t)= (J*xi+(delta/2)*((((2*B_grid)-J)*Ds)+J))*Pr;     % Gives the random sample (row) for each trajectory for each input parameter (column)   %  
    Ps_A(:,:,t)=Ps(:,1:nnz(A),t);                               % The trajectory is divided over two matrices: one for matrix A and one for matrix B %
    Ps_B(:,:,t)=Ps(:,nnz(A)+1:k,t);
end

R=0.1;                  % Range=10%%
A_upper=zeros(size(A)); % Scale each uncertain parameter to [0,1]. Assume that the uncerainty ranges are "10%", e.g. 50 becomes: [45 - 55] %
B_upper=zeros(size(B)); % This is done by creating the upper (A_upper, B_upper) and lower (A_lower, B_lower) valued matrices defined by the uncertainy ranges. 
A_lower=zeros(size(A));
B_lower=zeros(size(B));

for i=1:size(A,1);      % The upper values of the A-matrix %
    for j=1:size(A,1);
        if A(i,j)>0;
            A_upper(i,j)=A(i,j)+(A(i,j)*R);   
        else
            A_upper(i,j)=A(i,j)-(A(i,j)*R);
        end
    end
end

for i=1:size(A,1);      % The lower values of the A-matrix %
    for j=1:size(A,1);
        if A(i,j)>0;
            A_lower(i,j)=A(i,j)-(A(i,j)*R);  
        else
            A_lower(i,j)=A(i,j)+(A(i,j)*R);
        end
    end
end

for i=1:size(B,1);      % These are the upper values of the B-matrix %
    for j=1:size(B,2);
        if B(i,j)>0;
            B_upper(i,j)=B(i,j)+(B(i,j)*R) ; 
        else
            B_upper(i,j)=B(i,j)-(B(i,j)*R);
        end
    end
end

for i=1:size(B,1);      % The lower values of the B-matrix %
    for j=1:size(B,2);
        if B(i,j)>0;
            B_lower(i,j)=B(i,j)-(B(i,j)*R) ; 
        else
            B_lower(i,j)=B(i,j)+(B(i,j)*R);
        end
    end
end

h_A=find(A);          % The non-zero elements of A are given by the linear index h_A (e.g. (1,1)=1 and (2,1)=2 etc. %
h_B=find(B);          % The non-zero elements of G are given by the linear index h_B %

A_store=zeros(size(A));   % Pre-allocate for speed, A_step is the new matrix through the hypercube, of the same size as A % 
A_step_norm=zeros(size(A));
A_step=zeros(size(A));

B_store=zeros(size(B));
B_step_norm=zeros(size(B));
B_step=zeros(size(B));

g=zeros(k1,z);

% Calculate LCA result g for each point in the hypercube of all trajectories %

for t=1:z;                 
    for j=1:k1;
        A_store(h_A)=Ps_A(j,:,t);                                           % The non-zero elements of A are replaced by the subsequent points of the trajectory for each trajectory %
        A_step_norm(:,:,j,t)=A_store;
        A_step(:,:,j,t)=(A_step_norm(:,:,j,t)).*(A_upper-A_lower)+A_lower;  % The new parameters in A_step are scaled to their original interval
        
        B_store(h_B)=Ps_B(j,:,t);                                           % The same is done for the B-matrix %
        B_step_norm(:,:,j,t)=B_store;
        B_step(:,:,j,t)=(B_step_norm(:,:,j,t)).*(B_upper-B_lower)+B_lower;
        
        g(j,t)=B_step(:,:,j,t)*(A_step(:,:,j,t)\f);                         % Calculate the total environmental impact for each step of the trajectory for each trajectory %
    end
end

        
% Elementary effect %
ind=zeros(t,k);
for t=1:z; 
    for i=1:k;
        [ind(t,i)]=find(diff(Ps(:,i,t)));       % find the index (ind) for which the trajectory makes a step in the grid for each input parameter %
    end
end

% Calculate the elementary effects % 

EE=zeros(z,k);

for t=1:z;  
    for i=1:k;
         if (Ps(ind(t,i),i,t) > Ps(ind(t,i)+1,i,t)); 
            EE(t,i)= (g(ind(t,i),t) - g(ind(t,i)+1,t))/(2/3); 
         else 
            EE(t,i)= (g(ind(t,i)+1,t) - g(ind(t,i),t))/(2/3);
         end
    end
end
 

% Calculate sensitivity measures %

mu= mean(EE);                   % Mean 
mu_star=mean(abs(EE));          % Absolute mean
sd=std(EE);                     % Standard deviation

[L,IL]=sort(mu_star);
[K,IK]=sort(mu_star, 'descend');
[S,IS]=sort(sd,'descend');
figure (1)
bar(L);
set(gca,'XTickLabel',{IL})      %Only the nonzero parameters are included
xlabel('parameter')
ylabel('Sensitivity index (\mu*)')



