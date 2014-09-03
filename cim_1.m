clc
close all
clear all

N=1400;
M=200;%the length between interface_1 & interface_2
epi_minus=1;
epi_plus=30;
mu_minus=1;
mu_plus=30;% wood 4 &1.00000043
alpha=1/2; %put interface at middle
beta=1/2;
pl_gap=50;
frequency_check=2;

h=2*pi/(N-1);
Dx=zeros(N);
O=zeros(N/2);
I=eye(N);
I_f=eye(N);
value_Ez=zeros(N,1);

tic;

rho_plus_1=epi_plus/(alpha*epi_plus+beta*epi_minus);
rho_minus_1=epi_minus/(alpha*epi_plus+beta*epi_minus);
rho_plus_2=epi_minus/(alpha*epi_minus+beta*epi_plus);
rho_minus_2=epi_plus/(alpha*epi_minus+beta*epi_plus);
for i=1:N-1
   Dx(i*N+i)=1;
end
Dx=Dx+transpose(Dx);
Dx=Dx-2*I;
Dx(N*(N-1)+1)=1;
Dx(N)=1;
Dx= -epi_minus*(1/h^2)*Dx;
%---assume interface_1 on the N/2~N/2+1
Dx(N/2,:)=0;
Dx(N/2+1,:)=0;
    % on the grid_N/2
    Dx(N/2,N/2-1)= -epi_minus/h^2;
    Dx(N/2,N/2)= epi_minus*(1+rho_plus_1)/h^2;
    Dx(N/2,N/2+1)= -epi_minus*rho_plus_1/h^2;
    %on the grid_N/2+1
    Dx(N/2+1,N/2)= -epi_plus*rho_minus_1/h^2;
    Dx(N/2+1,N/2+1)= epi_plus*(1+rho_minus_1)/h^2;
    Dx(N/2+1,N/2+2)= -epi_plus/h^2;
    
% medium inside permitivity
for v=N/+2:N/2+M-1 
    Dx(v,:)= (epi_minus/epi_plus)^(-1)*Dx(v,:);
end
 %---assume interface_2 on the N/2+M~N/2+M+1
 Dx(N/2+M,:)=0;
 Dx(N/2+M+1,:)=0;
     %on the grid N/2+M
     Dx(N/2+M,N/2+M-1)= -epi_plus/h^2;
     Dx(N/2+M,N/2+M)= epi_plus*(1+rho_plus_2)/h^2;
     Dx(N/2+M,N/2+M+1)= -epi_plus*rho_plus_2/h^2;
     %on the grid N/2+M+1
     Dx(N/2+M+1,N/2+M)= -epi_minus*rho_minus_2*h^2;
     Dx(N/2+M+1,N/2+M+1)= epi_minus*(1+rho_minus_2)/h^2;
     Dx(N/2+M+1,N/2+M+2)= -epi_minus/h^2;
     
 I_f=mu_minus*I_f ;
 
for  t=N/2+1:N/2+M
    I_f(t,t)=(mu_plus/mu_minus)*I_f(t,t);
end

%--eigenvalue
[vec_Ez,lambda_Ez]=eig(Dx,I_f);

for s=1:N
    value_Ez(s)=sqrt(lambda_Ez(s,s));
end
val_want_Ez=find(value_Ez <frequency_check & value_Ez > -frequency_check);
value_Ez_show=value_Ez(value_Ez <frequency_check & value_Ez > -frequency_check);

criterion=zeros(length(val_want_Ez),4);
criterion_divide=zeros(length(val_want_Ez),2);
k_cnt=1;
for i=1:length(val_want_Ez)
    figure;
%     subplot(4,2,1:2);
    plot(real(vec_Ez(:,val_want_Ez(k_cnt))));
    grid on;
%     subplot(3,2,3);
%     plot(real(vec_Ez(:,val_want_Ez(k_cnt))));
%     axis([N/2-pl_gap N/2 -1 1]);
%     grid on;
%     subplot(3,2,4);
%     plot(real(vec_Ez(:,val_want_Ez(k_cnt))));
%     axis([N/2 N/2+pl_gap -1 1]);
%     grid on;
%     subplot(3,2,5);
%     plot(real(vec_Ez(:,val_want_Ez(k_cnt))));
%     axis([N/2+M-pl_gap N/2+M -1 1]);
%     grid on;
%     subplot(3,2,6);
%     plot(real(vec_Ez(:,val_want_Ez(k_cnt))));
%     axis([N/2+M N/2+M+pl_gap -1 1]);
%     grid on;
    criterion(i,1)=epi_minus*(vec_Ez(N/2,val_want_Ez(k_cnt))-vec_Ez(N/2-1,val_want_Ez(k_cnt)));
    criterion(i,2)=epi_plus*(vec_Ez(N/2+1,val_want_Ez(k_cnt))-vec_Ez(N/2,val_want_Ez(k_cnt)));
    criterion(i,3)=epi_plus*(vec_Ez(N/2+M,val_want_Ez(k_cnt))-vec_Ez(N/2+M-1,val_want_Ez(k_cnt)));
    criterion(i,4)=epi_minus*(vec_Ez(N/2+M+1,val_want_Ez(k_cnt))-vec_Ez(N/2+M,val_want_Ez(k_cnt)));
    k_cnt=k_cnt+1;
end
toc;