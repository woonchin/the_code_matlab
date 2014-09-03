function [omega,s1,s2]=wave1D(a,time,nx)
% a= the length of interval
% time= the total time interval 
% nx=the number of subunterval in (o,a)

% omega= angular frequency
% s1=the complex Fourier transf. at x=a/5
% s2=the complex Fourier transf. at x=a/2

f0= randn(nx+1,1);
f0(1)=0;
f0(nx+1)=0;

f1=randn(nx+1,1);
f1(1)=0;
f1(nx+1)=0;

dx=a/nx;
d2tmax=1.9*dx;

ntime=round(time/d2tmax+1);
dt=time/(2*ntime);

A=spalloc(nx+1,nx+1,3*(nx+1));
for i=2:nx
    A(i,i)=2*(1-(dt/dx)^2);
    A(i,i+1)=(dt/dx)^2;
    A(i,i-1)=(dt/dx)^2;
end

for itime=1:ntime
    
    f0=A*f1-f0;
    sign1(2*itime-1)=f0(round(1+nx/5));
    sign2(2*itime-1)=f0(round(1+nx/2));
    
    f1= A*f0-f1;
    sign1(2*itime)=f1(round(1+nx/5));
    sign2(2*ntime)=f1(round(1+nx/2));
   
end

spec1=fft(sign1);
spec2=fft(sign2);

s1(1:ntime)=spec1(1:ntime);
s2(1:ntime)=spec2(1:ntime);

omega = (2*pi/time)*linspace(0,ntime-1,ntime);