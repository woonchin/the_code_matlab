% [a,id]=sort(abs(val_show));
% id
% i = id(find(a>.003,1,'first'))
% a(i)
a=0;
b=600;
for i=1:length(val_want)
s=val_want(i);
pause(.5);
subplot(4,1,1);
plot(real(vec_Hx(:,s)));
grid on;
axis([a b -1 1])
title('real-Ex');
subplot(4,1,2);
plot(real(vec_Hy(:,s)));
grid on;
axis([a b -1 1])
title('Ey');
subplot(4,1,3);
plot(real(vec_Ez(:,s)));
grid on;
axis([a b -1 1])
title('real-Hx');
subplot(4,1,4);
plot(real(vec_Hz(:,s)));
grid on;
axis([a b -1 1])
title('Hz');
end