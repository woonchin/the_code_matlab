function [val,err]=ddet0822(epi1,mu1,gamma1,epi2,mu2,gamma2,kz,kx,w)
R1= w*( gamma1^2-epi1*mu1);
R2= w*( gamma2^2-epi2*mu2);
k1=sqrt(w^2*(gamma1+sqrt(epi1*mu1))^2-kz^2)*sqrt(w^2*(gamma1-sqrt(epi1*mu1))^2-kz^2);
k2=sqrt(w^2*(gamma2+sqrt(epi2*mu2))^2-kz^2)*sqrt(w^2*(gamma2-sqrt(epi2*mu2))^2-kz^2);

A1=[-exp(i*k1*pi),-exp(-i*k1*pi),0,0,exp(i*k2*pi),exp(-i*k2*pi),0,0];
A2=[exp(2*i*kx*pi),exp(2*i*kx*pi),0,0,-exp(2*i*k2*pi),-exp(-2*i*k2*pi),0,0];
A3=[0,0,-exp(i*k1*pi),-exp(-i*k1*pi),0,0,exp(i*k2*pi),exp(-i*k2*pi)];
A4=[0,0,exp(-2*i*kx*pi),exp(-2*i*kx*pi),0,0,-exp(2*i*k2*pi),-exp(-2*i*k2*pi)];
A5=[-i*exp(i*k1*pi)*k1*gamma1/R1,i*exp(-i*k1*pi)*k1*gamma1/R1,-exp(i*k1*pi)*k1*mu1/R1,exp(-i*k1*pi)*k1*mu1/R1,i*exp(i*k2*pi)*k2*gamma2/R2,-i*exp(-i*k2*pi)*k2*gamma2/R2,exp(i*k2*pi)*k2*mu2/R2,-exp(-i*k2*pi)*k2*mu2/R2 ];
A6=[i*exp(2*i*kx*pi)*k1*gamma1/R1,-i*exp(2*i*kx*pi)*k1*gamma1/R1,exp(-2*i*kx*pi)*k1*mu1/R1,-exp(-2*i*kx*pi)*k1*mu1/R1,-i*exp(2*i*k2*pi)*k2*gamma2/R2,i*exp(-2*i*k2*pi)*k2*gamma2/R2,-exp(2*i*k2*pi)*k2*mu2/R2,exp(-2*i*k2*pi)*k2*mu2/R2];
A7=[exp(i*k1*pi)*k1*epi1/R1,-exp(-i*k1*pi)*k1*epi1/R1,-i*exp(i*k1*pi)*k1*gamma1/R1,i*exp(-i*k1*pi)*k1*gamma1/R1,-exp(i*k2*pi)*k2*epi2/R2,exp(-i*k2*pi)*k2*epi2/R2,i*exp(i*k2*pi)*k2*gamma2/R2,-i*exp(-i*k2*pi)*k2*gamma2/R2];
A8=[-exp(2*i*kx*pi)*k1*epi1/R1,exp(2*i*kx*pi)*k1*epi1/R1,i*exp(-2*i*kx*pi)*k1*gamma1/R1,-i*exp(-2*i*kx*pi)*k1*gamma1/R1,exp(2*i*k2*pi)*k2*epi2/R2,-exp(-2*i*k2*pi)*k2*epi2/R2,-i*exp(2*i*k2*pi)*k2*gamma2/R2,i*exp(-2*i*k2*pi)*k2*gamma2/R2];
A=[A1;A2;A3;A4;A5;A6;A7;A8];
val=det(A);
err=norm(val,2);
end