clc
close all
clear all

N=1800;
M=200;%the length between interface_1 & interface_2
epi_minus=1;
epi_plus=.06;
mu_minus=1;
mu_plus=.06;
alpha=1/2; %put interface at middle
beta=1/2;

h=2*pi/(N-1);
Dx=zeros(N);
O=zeros(N/2);
I=eye(N);
I_half=eye(N/2); % interface on the half
value_Ez=zeros(N,1);

tic;
%interface_1 corresponding coefficients
epi_head= (beta+beta^2)*(1/2+alpha)*epi_minus + (alpha+alpha^2)*(1/2+beta)*epi_plus;
rho_plus= epi_plus/epi_head;
rho_minus= epi_minus/epi_head;

a_i__1= (beta+beta^2)*rho_minus + alpha*(1+2*beta)*rho_plus;
a_i_0=  -(beta+beta^2)*rho_minus - (1+alpha)*(1+2*beta)*rho_plus;
a_i_1= (1+beta)^2*rho_plus;
a_i_2= -beta^2*rho_plus;

a_i1_1= (alpha+alpha^2)*rho_plus + beta*(1+2*alpha)*rho_minus;
a_i1_0= -(alpha+alpha^2)*rho_plus - (1+beta)*(1+2*alpha)*rho_minus;
a_i1__1= (1+alpha)^2*rho_minus;
a_i1__2= -alpha^2*rho_minus;

%interface_2 corresponding coefficients
epi_head_2=(beta+beta^2)*(1/2+alpha)*epi_plus + (alpha+alpha^2)*(1/2+beta)*epi_minus;
rho_plus_2= epi_minus/epi_head;
rho_minus_2= epi_plus/epi_head;

a_j__1= (beta+beta^2)*rho_minus_2 + alpha*(1+2*beta)*rho_plus_2;
a_j_0=  -(beta+beta^2)*rho_minus_2 - (1+alpha)*(1+2*beta)*rho_plus_2;
a_j_1= (1+beta)^2*rho_plus_2;
a_j_2= -beta^2*rho_plus_2;

a_j1_1= (alpha+alpha^2)*rho_plus_2 + beta*(1+2*alpha)*rho_minus_2;
a_j1_0= -(alpha+alpha^2)*rho_plus_2- (1+beta)*(1+2*alpha)*rho_minus_2;
a_j1__1= (1+alpha)^2*rho_minus_2;
a_j1__2= -alpha^2*rho_minus_2;

%---create 2 order derivative matrix
for i=1:N-1
   Dx(i*N+i)=1;
end
Dx=Dx+transpose(Dx);
Dx=Dx-2*I;
Dx(N*(N-1)+1)=1;
Dx(N)=1;
Dx= (-1/epi_minus)*(1/h^2)*Dx;
% medium inside permitivity
for v=N/+2:N/2+M-1 
    Dx(v,:)= (epi_minus/epi_plus)*Dx(v,:);
end

%---assume interface_1 on the N/2~N/2+1
Dx(N/2,:)=0;
Dx(N/2+1,:)=0;   
    % on the grid_N/2
    Dx(N/2,N/2-1)= -epi_minus*a_i__1/h^2;
    Dx(N/2,N/2)= -epi_minus*a_i_0/h^2;
    Dx(N/2,N/2+1)= -epi_minus*a_i_1/h^2;
    Dx(N/2,N/2+2)= -epi_minus*a_i_2/h^2;
    %on the grid_N/2+1
    Dx(N/2+1,N/2-1)= -epi_plus*a_i1__2/h^2;
    Dx(N/2+1,N/2)= -epi_plus*a_i1__1/h^2;
    Dx(N/2+1,N/2+1)=  -epi_plus*a_i1_0 /h^2;
    Dx(N/2+1,N/2+2)= -epi_plus*a_i1_1/h^2;
    

 %---assume interface_2 on the N/2+M~N/2+M+1
 Dx(N/2+M,:)=0;
 Dx(N/2+M+1,:)=0;
     % on the grid_N/2+M
    Dx(N/2+M,N/2+M-1)= -epi_plus*a_j__1/h^2;
    Dx(N/2+M,N/2+M)= -epi_plus*a_j_0/h^2;
    Dx(N/2+M,N/2+M+1)= -epi_plus*a_j_1/h^2;
    Dx(N/2+M,N/2+M+2)= -epi_plus*a_j_2/h^2;
    %on the grid_N/2+M+1
    Dx(N/2+1,N/2-1)= -epi_minus*a_j1__2/h^2;
    Dx(N/2+1,N/2)= -epi_minus*a_j1__1/h^2;
    Dx(N/2+1,N/2+1)=  -epi_minus*a_j1_0 /h^2;
    Dx(N/2+1,N/2+2)= -epi_minus*a_j1_1/h^2;
     
 I_f=[mu_minus*I_half O; O mu_minus*I_half];
 
for  t=N/2+1:N/2+M
    I_f(t,t)=(mu_plus/mu_minus)*I_f(t,t);
end

%--eigenvalue
[vec_Ez,lambda_Ez]=eig(Dx,I_f);

for s=1:N
    value_Ez(s)=sqrt(lambda_Ez(s,s));
end
val_want_Ez=find(value_Ez <2 & value_Ez > -2);
value_Ez_show=value_Ez(value_Ez <2 & value_Ez > -2);
toc;

k_cnt=1;
for i=1:length(val_want_Ez)
    figure;
    plot(real(vec_Ez(:,val_want_Ez(k_cnt))));
    grid on
    k_cnt=k_cnt+1;
end