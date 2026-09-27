function [u, omega] = VMD(signal, alpha, tau, K, DC, init, tol)
% VMD算法实现
% 输入参数：
%   signal: 输入信号
%   alpha: 带宽约束参数，控制模态的带宽
%   tau: 时间步长，影响收敛速度
%   K: 分解模态数
%   DC: 是否包含直流分量
%   init: 初始化中心频率，单位为Hz
%   tol: 收敛容差
% 输出参数：
%   u: 分解后的模态，每列代表一个模态
%   omega: 各模态的中心频率，单位为Hz

% 默认参数设置
if nargin < 7
    tol = 1e-5;  % 默认收敛容差
end
if nargin < 6
    init = 1;    % 默认初始化频率
end
if nargin < 5
    DC = 0;      % 默认不包含直流分量
end
if nargin < 4
    K = 2;       % 默认分解为2个模态
end
if nargin < 3
    tau = 0.1;   % 默认时间步长
end
if nargin < 2
    alpha = 1000; % 默认带宽约束参数
end

% 信号预处理
N = length(signal);  % 信号长度
PRI = 25e-3;        % 采样间隔，与main函数匹配
fs = 1/PRI;         % 采样频率，40Hz
t = (0:N-1)*PRI;    % 实际时间序列

% 频谱对称延拓
freqs = (0:N-1)*fs/N;  % 实际频率向量，范围[0, fs]
f_hat = fft(signal).';  % 信号的FFT

% 初始化中心频率（将Hz转换为归一化频率）
if length(init) == 1
    % 如果只提供一个初始频率，则均匀分布在[0,fs/2]范围内
    omega_plus = zeros(K,1);
    for i = 1:K
        omega_plus(i) = fs/2/K*(i-1);
    end
else
    % 如果提供多个初始频率，则直接使用（已经是Hz单位）
    omega_plus = init(:);
end

% 初始化其他变量
lambda_hat = zeros(N,1);  % 拉格朗日乘子
u_hat = zeros(N,K);       % 各模态的频域表示
u_hat_plus = zeros(N,K);  % 正频率部分
u_old_hat = zeros(N,K);   % 上一次迭代的模态

% 主循环
n = 0;                    % 迭代计数器
sum_uk = zeros(N,1);      % 除当前模态外其他模态的和
u_hat = repmat(f_hat,1,K)/K;  % 初始模态分配

% 预计算频率项
freq_sq = freqs.^2;       % 频率平方项

% 迭代更新
while 1
    n = n + 1;
    
    % 更新u_hat
    for k = K:-1:1
        % 计算除当前模态外其他模态的和
        sum_uk = sum(u_hat,2) - u_hat(:,k);
        
        % 更新u_hat
        denominator = 1 + alpha*(freq_sq - 2*freqs*omega_plus(k) + omega_plus(k)^2);
        u_hat(:,k) = (f_hat - sum_uk - lambda_hat/2)./(denominator.');
        
        % 更新omega
        u_hat_plus(:,k) = u_hat(:,k);
        
        % 计算中心频率（只考虑正频率部分）
        abs_u_plus = abs(u_hat_plus(:,k)).^2;
        % 使用正频率部分计算中心频率
        omega_plus(k) = sum(freqs(1:N/2).'.*abs_u_plus(1:N/2))/sum(abs_u_plus);
        
        % 确保中心频率在合理范围内
        omega_plus(k) = min(max(omega_plus(k), 0), fs/2);
    end
    
    % 更新拉格朗日乘子
    lambda_hat = lambda_hat + tau*(f_hat - sum(u_hat,2));
    
    % 判断收敛
    u_diff = norm(u_hat - u_old_hat,2)/norm(u_hat,2);
    if u_diff < tol || n > 500  % 添加最大迭代次数限制
        break;
    end
    u_old_hat = u_hat;
end

% 重构时域信号
u = zeros(N,K);
for k = 1:K
    u(:,k) = real(ifft(u_hat(:,k)));
end

% 输出中心频率（已经是Hz单位）
omega = omega_plus;
end 