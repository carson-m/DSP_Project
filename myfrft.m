function y = myfrft(x, a)
% x: time domain signal
% p: order of frft
% y: output

    a = mod(a,4); %将阶数a缩小到[0,4)范围
    N = length(x); %信号长度
    shft = rem((0:N-1)+fix(N/2),N)+1; %移位量，使a=+-1时的FrFT可用FFT计算
    sqrtN = sqrt(N); %归一化系数
    
    if(a == 0)
        y = x; %做0次DFT则序列不变
        return;
    end
    if(a == 1)
        y(shft,1) = fft(x(shft)) / sqrtN;
        return;
    end
    if(a == 2)
        y = flipud(x); %做两次DFT则序列反转
        return;
    end
    if(a == 3)
        y(shft,1) = ifft(x(shft)) * sqrtN;
        return;
    end

    if(a > 2) %缩小阶数a范围至[0,2)
        a = a - 2;
        x = flipud(x); %将高阶拆成2阶和a-2阶级联
    end

    %缩小阶数a范围至[0.5,1.5]
    if(a > 1.5)
        a = a - 1;
        x(shft,1) = fft(x(shft)) / sqrtN; %将高阶拆成1阶与a-1阶级联
    end
    if(a < 0.5)
        a = a + 1;
        x(shft,1) = ifft(x(shft)) * sqrtN; %将低阶拆成-1阶与a+1阶级联
    end

    x = [zeros(N-1,1) ; sincInterp(x) ; zeros(N-1,1)];

    tana2 = tan(a*pi/4);
    sina = sin(a*pi/2);
    % chirp premultiplication
    chrp = exp(-1i*pi/N*tana2/4.*(-2*N+2:2*N-2)'.^2);
    x = chrp.*x;
    
    % chirp convolution
    c = pi/N/sina/4;
    y = fconv(exp(1i*c*(-(4*N-4):4*N-4)'.^2),x);
    y = y(4*N-3:8*N-7)*sqrt(c/pi);
    
    % chirp post multiplication
    y = chrp.*y;
    
    % normalizing constant
    y = exp(-1i*(1-a)*pi/4)*y(N:2:end-N+1);
end