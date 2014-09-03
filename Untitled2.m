N=10;
h=figure;
for t=1:N
    subplot(2,1,1)
    plot(vec_Ey(:,t))
    axis tight
    axis([-1,N/5,-2,2]);
    subplot(2,1,2)
    plot(vec_Ez(:,t))
    axis tight
    axis([-1,N/5,-2,2]);
    Q(t)=getframe(h);
  
end
movie2avi(Q, 'myPeaks.avi','fps',8);

