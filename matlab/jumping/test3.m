

zz = [xx';xx'+1;nan(1,length(xx))];
t = [t;t;t];

zz = zz(:);
t = t(:);

figure;
tic
line(t,zz,'LineWidth',1)
toc