close all
t=1:0.01:10;
f=figure('Name',['Fig ' num2str(1)],'Visible','off');
for w = 1:10
    subplot(10,1,w);
    plot(t,sin(1.2^w*t));
end
set(f, 'Visible', 'on')
set(f, 'PaperUnits', 'inches');
set(f, 'PaperPosition', [0 0 8.5 11]);
saveas(f,'testfig1.pdf')
open('testfig1.pdf')