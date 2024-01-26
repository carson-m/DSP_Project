# DSP大作业

## Task1

设计程序`my_chirp`生成chirp信号

```matlab
function x = my_chirp(T, f0, f1, fs)
    B = f1 - f0; %带宽
    mu = B/T; %调频率
    t = 0 : 1/fs : T; %时间轴
    x = exp(1i*pi*(2*f0*t+mu*t.^2)); %输出信号
end
```

针对STFT分析设计MATLAB程序如下

```matlab
fs = 5e7;
t_sample = 0 : 1/fs : 1e-5;
x = my_chirp(1e-5,15e6,35e6,fs);
stft(x, fs,"Window",hann(42,'periodic'),"FrequencyRange","twosided");
```

得到如下时频联合分布

<img src=".\img\task1.png" alt="task1" style="zoom:50%;" />

可以看出该信号在0\~10us内从15MHz线性调频变为35MHz，符合要求。

## Task2

从频域分辨率角度，对于频域分布恒定的信号，想要将多个信号分辨开来，窗长越长越好。从时域分辨率的角度，对于频率分布随时间变化的信号，窗长越短，频率变化越小，信号就越接近平稳，就越能准确地描述各时间点的频谱特性。结合这两个角度，STFT的窗长要折中选择，不能太长也不能太短。若太长，则不同时间的频谱互相干扰，若太短，时域分辨率差。若太短，则对每个时间段内的频域分辨率不足。

<img src=".\img\task2.png" alt="task2" style="zoom:50%;" />

MATLAB的STFT默认使用Hanning窗，仅调整窗长，发现对本题参数与信号，当窗长为40~42左右可以取得较好的效果，在多数时间下可以检测出3个频率分量，且可以看出wave_data的频率变化特征。

## Task3

写本题的过程不太顺利，我从一开始就想把第一和第二小问一起写。

### 直接用STFT

我意识到STFT的时频分辨率不够，而且如果想直接利用STFT得到结果需要在二维平面上进行搜索，实现起来较为复杂。如果SNR比较高还容易些，如果SNR较低则会很容易误判或出现较大误差。

### Radon Wigner Transform

经过文献调研，一开始我找到了RWT(Radon Wigner Transform)方法。这种方法利用了LFM信号的WVD为冲激线谱，即
$$
x(t)=e^{j(\omega_0 t+0.5mt^2)}\\
W_x(t,\omega)=\int_{-\infty}^{\infty}x(t+\tau/2)x^*(t-\tau/2)e^{-j\omega\tau}d\tau=\delta(\omega-(\omega_0+mt))
$$
的特点，通过求原信号WVD的起止时间和起止频率得到所需参数。这种方法可以得到较为准确的结果，但需要在Radon变换的二维平面上搜索极值点，对斜率的求取精度要求较高，较为复杂。且WVD结果往往会有交叉项，使得对起止时间的判定更加困难。所以最终没有采用这种方法。

### Radon Ambiguity Transform

后来我查到了一种使用RAT(Radon Ambiguity Transform)的方法。这种方法利用了LFM信号的模棱函数(Ambiguity Function)为过原点线段的特点，将原来WVD的二维搜索问题转化为了一维搜索问题。只需关注Radon变换后过原点直线的斜率就好，极大简化了计算。(下图为某次实验中计算出的Ambiguity Function，可以看出有两个信号，根据其斜率可以确定原信号调频率$f_m$)

<img src=".\task3_2\untitled3.png" alt="untitled3" style="zoom:50%;" />

**模糊函数：**
$$
\chi (\tau ,f)=\int _{-\infty }^{\infty}s(t)s^{*}(t-\tau )e^{i2\pi ft}\,dt
$$
**确定斜率的方法：Radon变换**

计算图像中沿$\arg_{法线}=\theta$，距中心点距离为$d$的直线的积分，将图像中的每条不同的直线转化为变换域中不同的点。

具体地，我提取了Radon变换$d=0$的那一行，然后根据已知的待估计信号个数$N$取前$N$个极大值对应的$\theta$作为检测结果。真实的$\arg_{直线}=\phi=90\degree-\theta$。



用Radon变换求直线的角度为$\phi$，根据公式
$$
k=\tan\phi\\
f_m = B/T\\
k=B/f_s(归一化调频率)\\
$$
得到
$$
f_m=f_s^2*\tan\phi
$$
由此可以计算出调频率$f_m$。

一个细节问题：实际上由于MATLAB pambgfun函数输出频率间隔和$f_s$不同，代码中并非乘$f_s^2$，但原理上是一样的。

但这样也带来了问题，就是**信号的起止时间和起止频率无法得知**，仅能提取调频率$f_m$的大小。

虽然单用RAT无法得到要求，但可以用RAT得到的调频率在STFT的结果中简化搜索，使得基于STFT的直线搜索方法更加鲁棒。下面我是这样做的。

```matlab
k = fm./fs./(f(2)-f(1));
```

我通过这行代码计算出信号在STFT图像中对应的线段的预期斜率。有了这个，就可以**使用一条斜率固定的直线从上到下将STFT图像扫描一遍**，假设信号的调频率可以区分，则可以根据调频率唯一确定我们想要寻找的那条线段，从而进一步找到信号的起始时间和持续时间。

具体地，我先将STFT结果二值化，这样可以在一定程度上滤除掉一部分噪声的影响，得到图像如下。

<img src=".\task3_2\untitled5.png" alt="untitled5" style="zoom:50%;" />

接下来，根据RAT求得的斜率，用一条直线扫描这张图象，**每扫到一个位置记录下扫描到的点数，取点数最大的位置作为结果**。得到了这条线段所在的直线位置，就可以进一步在直线上搜索线段的起止点，进而求出参数。就如下图所示，黄色的细线为最终确定的直线位置。

<img src=".\task3_2\untitled6.png" alt="untitled6" style="zoom:50%;" />

但还有一个问题，那就是这样扫描会将没滤掉的噪声也计算进去。于是我对扫到的点的x坐标按顺序组成的向量做差分，这样就求出了每两个点之间的间隔大小。如果相邻两个扫到的点x坐标差别较大，则认为这两个点中至少有一个是噪声带来的。根据两点在整个向量中的位置可以决定该删掉哪一边。这样最终得到的结果就会和真实值比较接近了。现在可以计算出$t_0,T,f_m(调频率),B$，可以直接根据线段左端点得到起始频率$f_0$，这样也很容易实现。**这也是我最终采用的方法。**

### RAT-FrFT

我在调研时还发现了一种比较有意思的方法：分数阶傅里叶变换(FrFT)。

通过RAT可以确定FrFT的阶数a
$$
a = \frac{\cot^{-1}(-f_m/f_s^2)}{\frac{\pi}{2}}
$$
在此阶数下对原信号做FrFT，在理想情况下会得到一个近似冲激函数的图像，结果在u处极大，其他处很小。而这个u与待估计的LFM信号的中心频率$f_h$存在关系：
$$
f_h=|u*\csc(a*\pi/2)|
$$
这样就可以用FrFT得到原信号的中心频率$f_h$，再用
$$
f_0=f_h-\frac{f_mT}{2}
$$
即可计算出$f_0$。但在实际使用中我发现这样计算误差较大，所以最终没有采用这种方法。因为FrFT的结果受点数和阶数影响较大，如果求的阶数不够准确，或者点数不够，则频域可能不是一个理想的冲激函数。

## Task3.1结果

总时长26e(-6)s

采用一固定信号，真值：
$$
f_0 = 2e6\\
B = 1.8e7\\
T = 10e(-6)\\
t_0 = 1e(-6)
$$
各信噪比下蒙特卡洛实验所得RMSE

<img src=".\task3_1\t0_RMSE.png" alt="t0_RMSE" style="zoom:50%;" />

<img src=".\task3_1\T_RMSE.png" alt="T_RMSE" style="zoom:50%;" />

<img src=".\task3_1\f0_RMSE.png" alt="f0_RMSE" style="zoom:50%;" />

<img src=".\task3_1\B_RMSE.png" alt="B_RMSE" style="zoom:50%;" />

## Task3.2结果

已知两个LFM信号，真值为
$$
f_0 = [2e6,9e6]\\
B = [1.8e7,6e6]\\
T = [10e-6,11e-6]\\
t0 = [1e-6,8e-6]
$$
各信噪比下蒙特卡洛实验得RMSE

<img src=".\task3_2\t0_RMSE.png" alt="t0_RMSE" style="zoom:50%;" />

<img src=".\task3_2\T_RMSE.png" alt="T_RMSE" style="zoom:50%;" />

<img src=".\task3_2\f0_RMSE.png" alt="f0_RMSE" style="zoom:50%;" />

<img src=".\task3_2\B_RMSE.png" alt="B_RMSE" style="zoom:50%;" />