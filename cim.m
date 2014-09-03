clc
close all
clear all

N=400;
h=2*pi/(N-1);
Dx=zeros(N);
O=zeros(N/2);
I=eye(N);
I_half=eye(N/2);% interface on the half
value_Ez=zeros(N,1);

tic;

epi_minus=.5;
epi_plus=.5+.01;
mu_minus=.5;
mu_plus=.5+.01;
alpha=1/2; %put interface at middle
beta=1/2;
rho_plus=epi_plus/(alpha*epi_plus+beta*epi_minus);
rho_minus=epi_minus/(alpha*epi_plus+beta*epi_minus);

for i=1:N-1
   Dx(i*N+i)=1;
end
Dx=Dx+transpose(Dx);
Dx=Dx-2*I;
Dx(N*(N-1)+1)=1;
Dx(N)=1;
Dx=(1/h^2)*Dx;
for s=1:N/2-1
    Dx=(1/epi_minus)*Dx;
end
for t=N/2+2:N
    Dx=(1/epi_plus)*Dx;
end
%---assume interface on the N/2~N/2+1
Dx(N/2,:)=0;
Dx(N/2+1,:)=0;
% on the grid_N/2
Dx(N/2,N/2-1)=-epi_minus/h^2;
Dx(N/2,N/2)=epi_minus*(1+rho_plus)/h^2;
Dx(N/2,N/2+1)=-epi_minus*rho_plus/h^2;
%on the grid_N/2+1
Dx(N/2+1,N/2)=-epi_plus*rho_minus/h^2;
Dx(N/2+1,N/2+1)=epi_plus*(1+rho_minus)/h^2;
Dx(N/2+1,N/2+2)=-epi_plus/h^2;

%--eigenvalue
I_f=[mu_minus*I_half O; O mu_plus*I_half];
[vec_Ez,lambda_Ez]=eig(Dx,-I_f);

for s=1:N
    value_Ez(s)=sqrt((-1)*lambda_Ez(s,s));
end
val_want_Ez=find(value_Ez <2 & value_Ez > -2);
val_want_Ez_val=value_Ez(value_Ez <2 & value_Ez > -2);
toc;