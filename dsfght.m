
Y=[999 2999];
t=1;
while t<3
subplot(5,1,1);
plot(vec_Ey(:,Y(t)));
subplot(5,1,2);
plot(vec_Ey(:,Y(t)+1));
subplot(5,1,3);
plot(vec_Ey(:,Y(t)+2));
subplot(5,1,4);
plot(vec_Ey(:,Y(t)+3));
subplot(5,1,5);
plot(vec_Ey(:,Y(t)+4));
t=t+1;
end