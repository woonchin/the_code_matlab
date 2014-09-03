clear all;
p=.05;
r=10;
i=5;

for j=0:10
a[j]=1/(1+(r-1).*(1-p*j).^i);
end