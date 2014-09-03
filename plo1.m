cc=hsv(12);
close all
clc
for s=1:length(val_want)
    pause(1);
    plot(abs(criterion_err2{s}),'color',cc(s,:))
    hold on
end