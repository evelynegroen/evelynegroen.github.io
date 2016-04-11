%Quasi Monte Carlo sampling using matrix-based LCA % 
%Author: Evelyne Groen {evelyne [dot] groen [at] gmail [dot] com}

A=[10 0; -2 100];       %A-matrix
B=[1 10];               %B-matrix
f=[1000; 0];            %Scaling vector f

g_LCA=B*(A\f);          %Ceterministic answer (kg CO_2 emissions)

l=nnz(A)+nnz(B);        %Total number of nonzero parameters
h=find(A);              %nonzero elements of A
h2=find(B);             %nonzero elements of B

N=2^10;                  %Sample size. NB: should be a power of 2 in case of QMCS for uniform distributed samples
CV=0.05;                 %CV: coefficient of variation: set to 5%

pAm=zeros(size(A));
pAm2=zeros(2,2,N);
pBm=zeros(size(B));
pBm2=zeros(1,2,N);

p=sobolset(l,'Skip', 1);
p=scramble(p,'MatousekAffineOwen');     %Scramble: otherwise each run is the same

pA=p(1:N,1:nnz(A));
pB=p(1:N,nnz(A)+1:l);

for i=1:N
    pAm(h)=pA(i,:);
    pAm2(:,:,i)=pAm;        %Quasi Monte Carlo realizations in A-matrix
end

for i=1:N
    pBm(h2)=pB(i,:);
    pBm2(:,:,i)=pBm;        %Quasi Monte Carlo realizations in B-matrix
end

for i=1:N; 
        for a=1:2; 
            for b=1:2;
                if abs(A(a,b))>0;
                    A2(a,b,i)=norminv(pAm2(a,b,i),A(a,b), abs(CV*A(a,b)));
                end
                if abs(B(b))>0;
                    B2(:,b,i)=norminv(pBm2(:,b,i),B(:,b), abs(CV*B(:,b)));
                end
            end
        end
        g(i)=B2(:,:,i)*(A2(:,:,i)\f);
end

meang=mean(g);
stdg=std(g);

figure(1)
hist(g)
hTitle=title('Quasi Monte Carlo sampling');
hXLabel=xlabel('kg CO_{2}e');
hYLabel=ylabel('Number of realizations');