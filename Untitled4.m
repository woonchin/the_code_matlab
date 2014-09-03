close all
clear all
clc

N=1500;
h=2*pi/(N-1);

Dif=zeros(N);
i=sqrt(-1);
O=zeros(N);
b=eye(N);
value_EzHy=zeros(N,1);
value_EyHz=zeros(N,1);

tic;
%-----------create difference matrix---------------
for j=1:N-1
    Dif(j*N+j)=1;   
end
Dif_n=transpose(Dif);
Dif=(1/(2*h))*(Dif+(-1)*Dif_n);
Dif((N-2)*N+1)=-1/(2*h);
Dif(2*N)=1/(2*h);
D=[-Dif O;O -Dif];

B=i*[O b;b O];  
[vec_EzHy,lambda_EzHy]=eig(D,B);
[vec_EyHz,lambda_EyHz]=eig(-D,B);
%--------------eigenvalue ------------------
for s=1:2*N
    value_EzHy(s)=lambda_EzHy(s,s);
    value_EyHz(s)=lambda_EyHz(s,s);
end
M1=mat2cell(vec_EzHy,[N N]);
M2=mat2cell(vec_EyHz,[N N]);
vec_Ez=M1{1,1};
vec_Hy=M1{2,1};
vec_Ey=M2{1,1};
vec_Hz=M2{2,1};

val_want_EyHz=find(value_EyHz <2 & value_EyHz > -2);
val_want_EzHy=find(value_EzHy <2 & value_EzHy > -2);

[EyHz_cnt,n]=size(val_want_EyHz);
[EzHy_cnt,n]=size(val_want_EzHy);
 

toc;