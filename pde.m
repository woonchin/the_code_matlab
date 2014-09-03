close all;
clear all;

n = input( 'n_partition=' ); 
h = 1/n ;
x = linspace(-1,1,n+1);
c = func(x,n,h);
p=ones(n,1);
a=-c\p;


a1 = zero(n +2,1);
a1(1)=0;
a1(n +2) = 0;
a1(2:n +1) = a;

plot(a1);

% h=(x_i-x_i-1)/2;
% x_ii = x-h;







