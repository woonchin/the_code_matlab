[X,Y,Z]=meshgrid(-1:.5:1;-1:.5:1;-1:.5:1;)
(sqrt(X.^2+Y.^2)-2).^2+Z.^2-4=0;
mesh(X,Y,Z);