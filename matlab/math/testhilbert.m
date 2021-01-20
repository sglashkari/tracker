clc; 
close all;  
fs = 30e3;
t = 0:1/fs:1; 

x = cos(2*pi*5*t) + sin(2*pi*8*t) + cos(2*pi*10*t);

subplot(3,1,1)
plot(t,x)
xlim([0 1])

y = hilbert(x);
subplot(3,1,2)
plot(t,real(y),t,imag(y))
xlim([0 1])
legend('real','imaginary')
title('hilbert Function')

phi = angle(-y);
subplot(3,1,3)
plot(t,phi)
xlim([0 1])

