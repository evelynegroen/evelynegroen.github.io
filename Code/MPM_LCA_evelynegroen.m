%Multiplier method using matrix-based LCA % 
%Author: Evelyne Groen {evelyne [dot] groen [at] gmail [dot] com}

%This code is created by: Reinout Heijungs and Sangwon Suh, source: "The computational
%structure of life cycle assessment", ISBN: 978-90-481-6041-9 (Print) 978-94-015-9900-9 (Online)

A=[10 0; -2 100];   %A-matrix
B=[1 10];           %B-matrix
f=[1000; 0];        %Functional unit vector f
g_LCA=B*(A\f);      %Deterministic answer

s=A\f;              %Scaling vector s

Lambda=B*inv(A);
dgdA=zeros(size(A,1));
dgdB=zeros(1,size(A,1));

for k=1:size(g_LCA), 
     for i=1:size (A, 1),
         for j=1:size(A,2),
             dgdA(i,j,k)=-Lambda(k,i)*s(j);
         end 
     end
 end; 
 
 for k=1:size(g_LCA);
     R(:,:,k)=A/g_LCA(k);
 end
 
for k=1:size(g_LCA);
    GammaA(:,:,k)=R(:,:,k).*dgdA(:,:,k);
end
 
 for k=1: size(g_LCA), 
     for i=1: size(B,1), 
         for j=1: size(B,2), 
             if i==k
                 dgdB(i,j,k)=s(j);
             else
                 dgdB(i,j,k)=0;
             end
         end
     end
 end; 

 for k=1:size(g_LCA);
     P(:,:,k)=B/g_LCA(k);
 end
 
for k=1:size(g_LCA);
    GammaB(:,:,k)=P(:,:,k).*dgdB(:,:,k);
end

