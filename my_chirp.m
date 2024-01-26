function x = my_chirp(T, f0, f1, fs)
    B = f1 - f0; %带宽
    mu = B/T; %调频率
    t = 0 : 1/fs : T; %时间轴
    x = exp(1i*pi*(2*f0*t+mu*t.^2)); %输出信号
end

