fig = figure('position',[100 100 850 600]);
x = 1:10;
for i = 1:10
subplot(2,1,1);
plot(x);
x = fliplr(x);
subplot(2,1,2);
plot(rand(1,10));
f(i) = getframe(fig);
end
close all
[h, w, p] = size(f(1).cdata);  % use 1st frame to get dimensions
hf = figure; 
% resize figure based on frame's w x h, and place at (150, 150)
set(hf, 'position', [150 150 w h]);
axis off
movie(hf,f);
mplay(f)