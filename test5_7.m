clc
clear all
close all

N=100;
h=2*pi/(N-1);

Dif=zeros(N);
d=eye(N); % for the Dif
I=eye(N);
O=zeros(N);
value_Ez=zeros(N,1);
tic;

%----------create difference matrix  backard-------
for i=1:N-1
   Dif(i*N+i)=1;
end
Dif=Dif+transpose(Dif);
Dif=Dif-2*I;
Dif(N*(N-1)+1)=1;
Dif(N)=1;
Dif=(1/h^2)*Dif;


% D=[Dif O;O Dif];
% B=-1*[I O;O I];  

%--------------eigenvalue ------------------
[vec_Ez,lambda_Ez]=eig(Dif,-I);
for s=1:N
    value_Ez(s)=sqrt(lambda_Ez(s,s));
end

val_want_Ez=find(value_Ez <2 & value_Ez > -2);
val_want_Ez_1=value_Ez(value_Ez<2 & value_Ez>-2);
toc;