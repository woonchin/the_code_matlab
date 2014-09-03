
for i=1:length(val_want)
s=val_want(i);
pause(1);
plot(real(vec_Ex(:,s)));
grid on
end