%% 通感一体化系统仿真主程序
clear all; close all; clc;

%% 系统参数设置
% 通信参数
N = 64;                % OFDM子载波数
CP = 16;               % 循环前缀长度
mod_order = 4;         % QPSK调制阶数
code_rate = 1/2;       % 信道编码率
num_symbols = 100;     % OFDM符号数量

% FMCW参数
B = 3.99e9;            % 带宽
T = 57e-6;             % 脉宽
K = B/T;               % 调频斜率
PRI = 25e-3;           % Chirp重复间隔
fs = 4e6;              % 采样率(相对于中频信号)
fc = 77e9;             % 载波频率
Ts = 1/fs;             % 采样周期
c = 3e8;               % 光速
lambda = c/fc;         % 波长
N_range = 512;         % 距离采样点数
N_mc = 1024;           % 脉冲数

% 人体参数
R = 2;                  % 目标距离(m)
breath_rate = 0.333;    % 呼吸频率(Hz)
heart_rate = 1.217;     % 心跳频率(Hz)

%% 通信信号生成与处理
% 生成随机比特流
msg = randi([0 1], N*log2(mod_order)/2*num_symbols, 1);

% 信道编码（卷积码）
trellis = poly2trellis(7, [171 133]);  % 卷积码生成多项式
encoded_msg = convenc(msg, trellis);  % 卷积编码

% 将二进制信息转换为QPSK符号
msg_reshape = reshape(encoded_msg, 2, [])';    % 每2位一组
msg_decimal = bi2de(msg_reshape, 'left-msb');  % 转换为十进制

% 将数据分配到各个OFDM符号
msg_decimal = reshape(msg_decimal, N, num_symbols);

% QPSK调制
modulated_data = pskmod(msg_decimal, mod_order, pi/4);

% OFDM调制
ofdm_symbols = zeros(N+CP, num_symbols);
for i = 1:num_symbols
    time_domain = ifft(modulated_data(:, i), N);
    ofdm_symbols(:, i) = [time_domain(end-CP+1:end); time_domain];
end

% 将OFDM符号串行化
ofdm_serial = reshape(ofdm_symbols, [], 1);

%% FMCW信号生成
% 生成FMCW发射信号
t_slow = (0:N_mc-1)*PRI;
t_slow = ones(N_range,1)*t_slow;
t_fast = 2*R/c+(-N_range/2:N_range/2-1)*Ts;
t_fast = t_fast.'*ones(1,N_mc);
t = t_slow+t_fast;

% 生成呼吸和心跳引起的相位调制
breath_sig = 0.012*sin(2*pi*breath_rate*t);  % 呼吸相位调制
heart_sig = 0.0005*sin(2*pi*heart_rate*t);   % 心跳相位调制

% 生成回波信号（考虑时延和相位变化）
tau = 2*(R+breath_sig+heart_sig)/c;  % 往返时延

% 直接仿真中频回波FMCW信号
fmcw = rectpuls(t_fast-tau-T/2,T).*exp(1j*2*pi*K*tau.*t_fast).*exp(1j*2*pi*fc*tau);%相干积累后回波矩阵
fmcw_serial = reshape(fmcw,[],1);

% 分时复用发射中频信号
tx=[ofdm_serial;fmcw_serial]; %暂不考虑复杂帧结构

%% 信道模拟
% 添加高斯白噪声
SNR = 5;  % 信噪比(dB)
rx = [awgn(tx(1:(N+CP)*num_symbols),SNR,"measured");awgn(tx((N+CP)*num_symbols+1:end),SNR,"measured")]; % 通信与感知信号功率不同

%% 通信信号处理
% 将接收到的OFDM符号并行化
noisy_ofdm = rx(1:(N+CP)*num_symbols);
rx_ofdm = reshape(noisy_ofdm, N+CP, num_symbols);

% OFDM解调
rx_freq = zeros(N, num_symbols);
for i = 1:num_symbols
    rx_ofdm_symbol = rx_ofdm(CP+1:end, i);
    rx_freq(:, i) = fft(rx_ofdm_symbol, N);
end

% QPSK解调
demod_data = pskdemod(rx_freq, mod_order, pi/4);

% 将QPSK符号转换为二进制信息
demod_data = de2bi(demod_data, 2, 'left-msb');  % 转换为二进制
demod_data = reshape(demod_data.', [], 1);    % 串行化

%信道解码
decoded_msg = vitdec(demod_data, trellis, 5, 'trunc', 'hard');

% 计算误码率
[~, ber] = biterr(msg, decoded_msg);
fprintf('误码率: %.2f%%\n', ber*100);

%% FMCW信号处理
noisy_fmcw = rx((N+CP)*num_symbols+1:end);
rx_fmcw = reshape(noisy_fmcw,N_range,N_mc);

% 距离维FFT
range_fft = fftshift(fft(fftshift(rx_fmcw)));
r = (-fs/2:fs/N_range:fs/2-fs/N_range)/K*c/2;
figure
mesh(1:N_mc,r,abs(range_fft))
xlabel('Chirp number')
ylabel('Range(m)')

% 相位提取和解缠绕
[~,max_id] = max(abs(range_fft));
phase = angle(range_fft(max_id+N_range*(0:N_mc-1)));
figure
plot(phase)
grid on;
title('Phase before unwrapping');
xlabel('Time(s)');
ylabel('Phase');
phase_unwrap = unwrap(phase)/pi*180;
figure
plot(phase_unwrap)
grid on;
title('Phase after unwrapping');
xlabel('Time(s)');
ylabel('Phase');

% 带通滤波
% 呼吸信号滤波 (0.1-0.5Hz)
[b_breath, a_breath] = butter(4, [0.1 0.5]/(1/PRI/2), 'bandpass');
breath_signal = filter(b_breath, a_breath, phase_unwrap);

% 心跳信号滤波 (0.8-2Hz)
[b_heart, a_heart] = butter(4, [0.8 2]/(1/PRI/2), 'bandpass');
heart_signal = filter(b_heart, a_heart, phase_unwrap);

% 时域分析
figure;
subplot(2,1,1);
plot((0:N_mc-1)*PRI, breath_signal);
grid on;
title('Time domain diagram of breath signal');
xlabel('Time(s)');
ylabel('Amplitude');

subplot(2,1,2);
plot((0:N_mc-1)*PRI, heart_signal);
grid on;
title('Time domain diagram of heartbeat signal');
xlabel('Time(s)');
ylabel('Amplitude');

% 频谱分析
[P_breath, f_breath] = pwelch(breath_signal, [], [], [], 1/PRI);
[P_heart, f_heart] = pwelch(heart_signal, [], [], [], 1/PRI);
figure
subplot(2,1,1);
plot(f_breath,abs(P_breath))
grid on;
title('Spectrum of breath signal');
xlabel('Doppler(Hz)')
ylabel('Amplitude')
subplot(2,1,2);
plot(f_heart,abs(P_heart))
grid on;
title('Spectrum of heartbeat signal');
xlabel('Doppler(Hz)')
ylabel('Amplitude')

%% 人体生命体征检测
%%%%%%%%%%%%%%%%%%  Method 1: 2D-FFT  %%%%%%%%%%%%%%%%%%
% 呼吸信号检测
breath_signal_fft = fftshift(fft(fftshift(breath_signal),[]));
% 心跳信号检测
heart_signal_fft = fftshift(fft(fftshift(heart_signal),[]));
doppler = (0:N_mc/2-1)/N_mc/PRI;
[~,bre_index] = max(abs(breath_signal_fft(N_mc/2+1:end)));
fre_breath_2dFFT = doppler(bre_index);
fprintf('Estimation of breath rate with 2D-FFT: %.3f Hz\n', fre_breath_2dFFT);
[~,heart_index] = max(abs(heart_signal_fft(N_mc/2+1:end)));
fre_heartbeat_2dFFT = doppler(heart_index);
fprintf('Estimation of heartbeat rate with 2D-FFT: %.3f Hz\n', fre_heartbeat_2dFFT);

figure
subplot(2,1,1);
plot(doppler,abs(breath_signal_fft(N_mc/2+1:end)))
grid on;
title('Spectrum of breath signal(2D-FFT)');
xlabel('Doppler(Hz)')
ylabel('Amplitude')
subplot(2,1,2);
plot(doppler,abs(heart_signal_fft(N_mc/2+1:end)))
grid on;
title('Spectrum of heartbeat signal(2D-FFT)');
xlabel('Doppler(Hz)')
ylabel('Amplitude')

%%%%%%%%%%%%%%%%%%  Method 2: MUSIC  %%%%%%%%%%%%%%%%%%
% 呼吸信号检测
[fre_breath_MUSIC, P_breath_MUSIC] = MUSIC(breath_signal, N_mc, 1/PRI);
fprintf('Estimation of breath rate with MUSIC: %.3f Hz\n', fre_breath_MUSIC);
% 心跳信号检测
[fre_heart_MUSIC, P_heart_MUSIC] = MUSIC(heart_signal, N_mc, 1/PRI);
fprintf('Estimation of heartbear rate with MUSIC: %.3f Hz\n', fre_heart_MUSIC);
figure
subplot(2,1,1);
plot(doppler,abs(P_breath_MUSIC(1:N_mc/2)))
grid on;
title('Spectrum of breath signal(MUSIC)');
xlabel('Doppler(Hz)')
ylabel('Amplitude')
subplot(2,1,2);
plot(doppler,abs(P_heart_MUSIC(1:N_mc/2)))
grid on;
title('Spectrum of heartbeat signal(MUSIC)');
xlabel('Doppler(Hz)')
ylabel('Amplitude')

%%%%%%%%%%%%%%%%%%  Method 3: VMD  %%%%%%%%%%%%%%%%%%
% 参数设置（用 vmd_ 前缀区分开前面 FMCW 部分已经用过的 K / tau，避免变量名混用）
alpha = 100;         % 带宽约束参数，降低以加快收敛
vmd_tau = 0.1;       % 时间步长，增加以加快收敛
vmd_K = 2;           % 分解模态数，设置为2以分离呼吸和心跳
DC = 0;              % 不包含直流分量
init = [0.3, 1.2];   % 初始化中心频率，根据呼吸和心跳频率设置
tol = 1e-5;          % 收敛容差，适当放宽以加快收敛

% 对带通滤波后的信号进行VMD分解
[u_breath, omega_breath] = VMD(breath_signal, alpha, vmd_tau, vmd_K, DC, init, tol);
[u_heart, omega_heart] = VMD(heart_signal, alpha, vmd_tau, vmd_K, DC, init, tol);

% 输出频率估计结果
fprintf('Estimation of breath rate with VMD: %.3f Hz\n', omega_breath(1));
fprintf('Estimation of heartbeat rate with VMD: %.3f Hz\n', omega_heart(2));

% 显示分解结果
figure;
% 呼吸信号分解
subplot(2,2,1);
plot((0:N_mc-1)*PRI, u_breath(:,1));
grid on;
title('Breath signal with VMD-Mode 1');
xlabel('Time(s)');
ylabel('Amplitude')

subplot(2,2,2);
plot((0:N_mc-1)*PRI, u_breath(:,2));
grid on;
title('Breath signal with VMD-Mode 2');
xlabel('Time(s)');
ylabel('Amplitude')

% 心跳信号分解
subplot(2,2,3);
plot((0:N_mc-1)*PRI, u_heart(:,1));
grid on;
title('Heartbeat signal with VMD-Mode 1');
xlabel('Time(s)');
ylabel('Amplitude')

subplot(2,2,4);
plot((0:N_mc-1)*PRI, u_heart(:,2));
grid on;
title('Heartbeat signal with VMD-Mode 2');
xlabel('Time(s)');
ylabel('Amplitude')

% 计算各模态的频谱
figure;
% 呼吸信号频谱
subplot(2,2,1);
[P1_breath, f1_breath] = pwelch(u_breath(:,1), [], [], [], 1/PRI);
plot(f1_breath, abs(P1_breath));
grid on;
title('Breath signal with VMD-Mode 1');
xlabel('Doppler(Hz)')
ylabel('Amplitude')

subplot(2,2,2);
[P2_breath, f2_breath] = pwelch(u_breath(:,2), [], [], [], 1/PRI);
plot(f2_breath, abs(P2_breath));
grid on;
title('Breath signal with VMD-Mode 2');
xlabel('Doppler(Hz)')
ylabel('Amplitude')

% 心跳信号频谱
subplot(2,2,3);
[P1_heart, f1_heart] = pwelch(u_heart(:,1), [], [], [], 1/PRI);
plot(f1_heart, abs(P1_heart));
grid on;
title('Heartbeat signal with VMD-Mode 1');
xlabel('Doppler(Hz)')
ylabel('Amplitude')

subplot(2,2,4);
[P2_heart, f2_heart] = pwelch(u_heart(:,2), [], [], [], 1/PRI);
plot(f2_heart, abs(P2_heart));
grid on;
title('Heartbeat signal with VMD-Mode 2');
xlabel('Doppler(Hz)')
ylabel('Amplitude')







