function [freq_est, P_music] = MUSIC(signal, N, fs)
% MUSIC算法实现
% 输入参数：
%   signal: 输入信号
%   N: 信号长度
%   fs: 采样频率
% 输出参数：
%   freq_est: 估计的频率
%   P_music: MUSIC谱

% 参数设置
p = 1;          % 信号源数量（呼吸和心跳）
M = 128;         % 协方差矩阵的阶数，M越大，估计越精确
N_fft = N;      % FFT点数

% 构建数据矩阵
L = N - M + 1;  % 快拍数
X = zeros(M, L);
for i = 1:L
    X(:, i) = signal(i:i+M-1);
end

% 计算协方差矩阵
R = X * X' / L;

% 特征值分解
[V, D] = eig(R);
[D, I] = sort(diag(D), 'descend');
V = V(:, I);

% 噪声子空间
Vn = V(:, p+1:end);

% 计算MUSIC谱
f = (0:N_fft-1)/N_fft * fs;
P_music = zeros(1, N_fft);
for i = 1:N_fft
    a = exp(-1j*2*pi*(0:M-1)'*f(i)/fs);
    P_music(i) = 1 / (a' * (Vn * Vn') * a);
end

% % 归一化
% P_music = abs(P_music) / max(abs(P_music));

% 寻找峰值
[~, locs] = max(abs(P_music(1:N/2)));
freq_est = f(locs(1:p));

end 