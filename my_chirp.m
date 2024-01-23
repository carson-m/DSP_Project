function x = my_chirp(T, f0, f1, fs)
    B = f1 - f0;
    mu = B/T;%调频率
    t = 0:1/fs:T;
    % n=round(T*fs);%采样点个数
    % t = linspace(0,T,n);
    x = exp(2i*pi*(f0*t+0.5*mu*t.^2));
end

