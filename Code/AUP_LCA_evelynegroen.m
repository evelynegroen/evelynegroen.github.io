%Analytical uncertainty propagation using matrix-based LCA % 
%Author: Evelyne Groen {evelyne [dot] groen [at] gmail [dot] com}
%This code is based on Heijungs and Suh (2002); http://dx.doi.org/10.1007/978-94-015-9900-9

A=[10 0; -2 100];       %A-matrix
B=[1 10];               %B-matrix
f=[1000; 0];            %Scaling vector f

CV=0.05;                %CV: coefficient of variation
varA=abs(A*CV).^2;      %Variance of parameters in A-matrix
varB=abs(B*CV).^2;      %Variance of parameters in B-matrix

s=A\f;
g=B*s;

tic
Lambda=B*inv(A);

dgdA=zeros(size(A,1));
dgdB=zeros(1,size(A,1));
vargA=zeros(1,size(A,1));
vargB=zeros(1,size(A,1));
varg=zeros(1);

for k=1:size(g), 
     for i=1:size (A, 1),
         for j=1:size(A,2),
             dgdA(i,j,k)=-Lambda(k,i)*s(j);
         end 
     end
 end; 
  
 
 for k=1: size(g), 
     for i=1: size(B,1), 
         for j=1: size (B,2), 
             if i==k
                 dgdB(i,j,k)=s(j);
             else
                 dgdB(i,j,k)=0;
             end
         end
     end
 end; 

for k=1:size(g)
    vargA(k)=sum(sum((dgdA(:,:,k).^2).*varA));
end
 
for k=1:size(g)
     vargB(k)=sum(sum((dgdB(:,:,k).^2).*varB));
 end
 

for k=1:size(g),
    varg(k)=vargA(k)+vargB(k);
end

stdg=sqrt(varg);
