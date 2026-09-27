% FMCW雷达呼吸心跳检测仿真
clear all; close all; clc;

% 参数设置
fc = 24e9;          % 载波频率 24GHz
B = 1e9;            % 带宽 1GHz
T = 0.1;            % 扫频周期 100ms
fs = 2e6;           % 采样率 2MHz
c = 3e8;            % 光速
K = B/T;            % 调频斜率
N_scan = 100;       % 扫描次数

% 目标参数
d0 = 1;             % 初始距离 1m
A_breath = 0.01;    % 呼吸幅度 1cm
A_heart = 0.001;    % 心跳幅度 1mm
f_breath = 0.2;     % 呼吸频率 0.2Hz
f_heart = 1.2;      % 心跳频率 1.2Hz

% 时间轴
t_scan = 0:1/fs:T-1/fs;    % 单个扫描周期的时间
t_total = 0:1/fs:N_scan*T-1/fs;  % 总时间
N = length(t_scan);

% 初始化存储数组
phase_history = zeros(1, N_scan);
breath_signal = zeros(1, N_scan);
heart_signal = zeros(1, N_scan);

% 对每个扫描周期进行处理
for i = 1:N_scan
    % 当前时间点
    t_current = t_scan + (i-1)*T;
    
    % 生成目标运动
    d_breath = A_breath*sin(2*pi*f_breath*t_current);
    d_heart = A_heart*sin(2*pi*f_heart*t_current);
    d = d0 + d_breath + d_heart;
    
    % 生成发射信号
    s_t = cos(2*pi*fc*t_scan + pi*K*t_scan.^2);
    
    % 生成接收信号
    tau = 2*d/c;
    s_r = cos(2*pi*fc*(t_scan-tau) + pi*K*(t_scan-tau).^2);
    
    % 混频
    s_if = s_t.*s_r;
    
    % 距离FFT
    N_fft = 2^nextpow2(N);
    f = (-N_fft/2:N_fft/2-1)*fs/N_fft;
    S_if = fftshift(fft(s_if, N_fft));
    
    % 找到目标距离
    [~, idx] = max(abs(S_if));
    
    % 提取相位
    phase_history(i) = angle(S_if(idx));
end

% 相位解缠绕
phase_unwrap = unwrap(phase_history);

% 带通滤波
% 呼吸信号滤波 (0.1-0.5Hz)
[b_breath, a_breath] = butter(4, [0.1 0.5]/(1/T/2), 'bandpass');
breath_signal = filter(b_breath, a_breath, phase_unwrap);

% 心跳信号滤波 (0.8-2Hz)
[b_heart, a_heart] = butter(4, [0.8 2]/(1/T/2), 'bandpass');
heart_signal = filter(b_heart, a_heart, phase_unwrap);

% 绘图
figure;
subplot(3,1,1);
plot((0:N_scan-1)*T, phase_unwrap);
title('解缠绕后的相位');
xlabel('时间(s)');
ylabel('相位(rad)');

subplot(3,1,2);
plot((0:N_scan-1)*T, breath_signal);
title('提取的呼吸信号');
xlabel('时间(s)');
ylabel('幅度');

subplot(3,1,3);
plot((0:N_scan-1)*T, heart_signal);
title('提取的心跳信号');
xlabel('时间(s)');
ylabel('幅度');

% 频谱分析
figure;
subplot(2,1,1);
[P_breath, f_breath] = pwelch(breath_signal, [], [], [], 1/T);
plot(f_breath, 10*log10(P_breath));
title('呼吸信号频谱');
xlabel('频率(Hz)');
ylabel('功率谱密度(dB/Hz)');

subplot(2,1,2);
[P_heart, f_heart] = pwelch(heart_signal, [], [], [], 1/T);
plot(f_heart, 10*log10(P_heart));
title('心跳信号频谱');
xlabel('频率(Hz)');
ylabel('功率谱密度(dB/Hz)'); 