%% 通感一体化系统仿真主程序
clear all; close all; clc;

%% 系统参数设置
% 通信参数
N = 128;                % OFDM子载波数
CP = 16;               % 循环前缀长度
mod_order = 4;         % QPSK调制阶数
code_rate = 1/2;       % 信道编码率
num_symbols = 1000;     % OFDM符号数量

% 噪声水平
SNR = 0:3:15;  % 信噪比(dB)
MaxIter = 2e2;

%% 无信道编码
% 生成随机比特流
msg = randi([0 1], N*log2(mod_order)*num_symbols, 1);

% 将二进制信息转换为QPSK符号
msg_reshape = reshape(msg, 2, [])';    % 每2位一组
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
ber1 = zeros(MaxIter,length(SNR));
for index=1:length(SNR)
    for m=1:MaxIter
        rx = awgn(ofdm_serial,SNR(index),"measured");

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
        demod_msg = de2bi(demod_data, 2, 'left-msb');  % 转换为二进制
        demod_msg = reshape(demod_msg.', [], 1);    % 串行化

        % 计算误码率
        [~, ber_temp] = biterr(msg, demod_msg);
        ber1(m,index)=ber_temp;
    end
end
ber1=mean(ber1);

%% 有信道编码
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
ber2 = zeros(MaxIter,length(SNR));
for index=1:length(SNR)
    for m=1:MaxIter
        rx = awgn(ofdm_serial,SNR(index),"measured");

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
        [~, ber_temp] = biterr(msg, decoded_msg);
        ber2(m,index)=ber_temp;
    end
end
ber2=mean(ber2);
ber2(ber2==0)=eps;

figure
semilogy(SNR,ber1)
hold on; grid on; box on;
semilogy(SNR,ber2)
legend('Without convolutional code','With convolutional code')
xlabel('SNR/dB')
ylabel('BER')







