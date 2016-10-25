%Procedure:     Uncertainty propagation for matrix-based LCA 
%Method:        Fuzzy interval arithmetic
%Author:        Evelyne Groen {evelyne [dot] groen [at] gmail [dot] com}
%Last update:   20/10/2016 
%Toolbox:       none

A=[10 0; -2 100];       %A-matrix
B=[1 10];               %B-matrix
f=[1000; 0];            %Functional unit vector f

g_LCA=B*(A\f);          %Determinstic answer (kg CO_2)

R=0.1;                  %R:range, set to 10%
A_lower=A-A*R;          %Lower bound of parameters in the A-matrix
A_upper=A+A*R;          %Upper bound

B_lower=B-B*R;          %Lower bound of parameters in the B-matrix
B_upper=B+B*R;          %Upper bound

g_upper=B_upper*(A_lower\f); % Boundaries of the output
g_mean=B*(A\f);
g_lower=B_lower*(A_upper\f);

a=11;                   % Number of alpha cuts
g_up=zeros(1,a);        % Pre-allocate for speed
g_low=zeros(1,a);

for k=1:a;
    A1(:,:) =((k-1)/10)*(A(:,:)-A_lower(:,:))+A_lower(:,:); %Lower A-matrix cuts
    A2(:,:)=(-(k-1)/10)*(A_upper(:,:)-A(:,:))+A_upper(:,:); %Upper A-matrix cuts
         
    B1(:,:)=(-(k-1)/10)*(B_upper(:,:)-B(:,:))+B_upper(:,:); %Upper B-matrix cuts
    B2(:,:)=((k-1)/10)*(B(:,:)-B_lower(:,:))+B_lower(:,:);  %Lower B-matrix cuts
   
    g_up(:,k)=B1(:,:)*(A1(:,:)\f);
    g_low(:,k)=B2(:,:)*(A2(:,:)\f); 
end

figure(1)
k=0:0.1:1; 
plot(g_up,k, '-o')
hold on
plot(g_low,k, '-o')
hTitle=title('Fuzzy interval arithmetic');
hXLabel=xlabel('kg CO_{2}e');
hYLabel=ylabel('Alpha-cut');