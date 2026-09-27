close all
clear
clc

% 呼吸心跳信号仿真
ts=1e-2;%相干积累后时间
N=1024;
t=0:ts:(N-1)*ts;
f=(0:N-1)/N/ts;
fb=0.333;
fh=1.1217;
x=sin(2*pi*fb*t)+sin(2*pi*fh*t);
y=(fft(fftshift(x)));
figure
subplot(2,1,1)
plot(t,x)
grid on;
title('Time domain diagram of mixed signal')
xlabel('Time(s)');
ylabel('Amplitude')
subplot(2,1,2)
plot(f(1:100),abs(y(1:100)))
grid on;
title('Spectrum of mixed signal')
xlabel('f(Hz)');
ylabel('Amplitude')






