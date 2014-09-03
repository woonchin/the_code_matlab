Z = peaks;
figure('Renderer','zbuffer');
subplot(1,2,1)
surf(Z);title('first plot')
axis tight;
set(gca,'NextPlot','replaceChildren');
subplot(1,2,2);
surf(-Z);title('second plot')
axis tight;
set(gca,'NextPlot','replaceChildren');
for jj = 1:20
    subplot(1,2,1);
    surf(sin(2*pi*jj/20)*Z,Z)
    subplot(1,2,2);
    surf( -sin(2*pi*jj/20)*Z,Z);
    F(jj) = getframe;
end
movie2avi(F, 'mymov.avi', 'Compression','none');