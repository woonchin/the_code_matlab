clear all
close all
clc
N=300;

a=2*pi;%unit cell
cnt=0;
ccnt_1=0;%08/14


alpha=1/2;%(x_head - x_i) /h assume put on the middle
beta=1/2; %(x_i+1-x_head)/h
epi_minus=1;
mu_minus=1;
gamma_minus=0;
kz=0;
epi_plus=20;
mu_plus=2;
gamma_plus=0.3;
check_frequency=1;
h=a/N;
%%  assume the interface_1 at N/2~N/2+1 and interface_2 at N/2+M ~N/2+M+1

tic;

Dif=zeros(N);
d=eye(N); % for the Dif
I=eye(N);
O_3n=zeros(3*N);
O=zeros(N);
eig_val=zeros(6*N,1);

%%  difference matrix A with Ghost state
for s=1:N-1
    Dif(s*N+s)=1;
end
D_E=(1/h)*(Dif-d);
D_E(N)=(1/h)*1;
D_H=(-1)*transpose(D_E);

%     D_E(N)=0;
%     D_H(N*(N-1)+1)=0;
%      M=N/2;
%     % interface_1&2
%     D_E(N/2,:)=0;
%     D_E(N/2+M,:)=0;
%     %for at N/2
%     D_E(N/2,N/2)=(1/h)*(-1/alpha);
%     D_E(N/2,N/2+1)=(1/h)*((1+beta)/alpha);
%     D_E(N/2,N/2+2)=(1/h)*(-beta/alpha);
%     %for at N/2+M
%     D_E(N/2+M,N/2+M)=(1/h)*(-1/alpha);
%     D_E(N/2+M,2)=(1/h)*((1+beta)/alpha);
%     D_E(N/2+M,3)=(1/h)*(-beta/alpha);
%      %interface_1&2
%     D_H(N/2+1,:)=0;
%     D_H(1,:)=0;
%     %for at N/2+1
%     D_H(N/2+1,N/2-1)=(1/h)*(alpha/beta);
%     D_H(N/2+1,N/2)=(1/h)*((-alpha-1)/beta);
%     D_H(N/2+1,N/2+1)=(1/h)*(1/beta);
%     %for at 2
%     D_H(1,N/2+M-1)=(1/h)*(alpha/beta);
%     D_H(1,N/2+M)=(1/h)*((-alpha-1)/beta);
%     D_H(1,1)=(1/h)*(1/beta);


for kx=0:.002:0.5
    D_E(N)=(1/h)*exp(i*kx*a); % bloch
    D_H(N*(N-1)+1)=-(1/h)*exp(-i*kx*a);
    
    
    A_1=[O -i*kz*I O; i*kz*I O -D_E; O D_E O];
    A_2=O_3n;
    A_3=O_3n;
    A_4=[O -i*kz*I O; i*kz*I O -D_H; O D_H O];
    A=[A_1 A_2; A_3 A_4];
    %% B
    B_1=gamma_minus*I;
    B_2=i*mu_minus*I;
    B_3=-i*epi_minus*I;
    B_4=gamma_minus*I;
    for s=N/2+1:N
        B_1(s,s)=gamma_plus*I(s,s);  %fix for gamma=0
        B_2(s,s)=i*mu_plus*I(s,s);
        B_3(s,s)=-i*epi_plus*I(s,s);
        B_4(s,s)=gamma_plus*I(s,s);
    end
    B_1=[B_1 O O;O B_1 O;O O B_1 ];
    B_2=[B_2 O O;O B_2 O;O O B_2 ];
    B_3=[B_3 O O;O B_3 O;O O B_3 ];
    B_4=[B_4 O O;O B_4 O;O O B_4 ];
    B=[B_1 B_2; B_3 B_4];
    
    %calculate eig
    [vec,lambda]=eig(A,B);
    for s=1:6*N
        eig_val(s)=lambda(s,s);
    end
    
    val_want=find(eig_val<=check_frequency & eig_val>0);
    val_show=eig_val(eig_val<=check_frequency & eig_val>0);
    %% error check 08/14
    for s=1:2:length(val_show) 
        [val_1,err_1]=ddet0822(epi_plus,mu_minus,gamma_minus,epi_plus,mu_plus,gamma_plus,kz,kx,val_show(s));
        error_1(ccnt_1+s,:)=[val_1,err_1];
    end
    ccnt_1=ccnt_1+length(val_show);
    %% plot setting
    Kx=kx*ones(length(val_show),1);
    for t=1:length(val_show)
        te(cnt+t,:)=[val_show(t),Kx(t)];
    end
    cnt=cnt+length(val_show);
end
toc
plot(real(te(:,2)),real(te(:,1)),'r.')
figure
plot(error_1(:,2))