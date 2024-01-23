%Program to Calculate CONTIONOUS FRACTIONAL FOURIER TRANSFORM
% By Ashish Soni,ATC Indore
tic
clc; 
clear all; 
close all;

T=input('enter the value of T - ');
t=-T:0.01:T;
X=input('enter input function - ');
disp('Program to Calculate Continous Fractional Fourier Transform')
U=input('enter the range of frequency - ');
a=input('Enter the rotation angle a (except 0) - ');

p=0;q=0; 
 for t=-T:0.01:T
  p=p+1;
    for u=-U:0.1:U;
        q=q+1;
%         x1(p,q)=(sqrt((1-j*cot(a))/(2*pi)))*exp(j*(((0.5*(u^2)+0.5*(t^2))*cot(a))-u*t*(1/sin(a))));
        x1(p,q)=(sqrt(1-j*cot(a*pi/2)))*exp(j*pi*((u^2)*cot(a*pi/2)+(t^2)*cot(a*pi/2)-2*u*t*(1/sin(a*pi/2))));
    end
    q=0;
 end
 
 
u=-U:0.1:U;
x2=0.01*(X*x1);

% x3=20*log10(abs(x2)./max(x2));
t=-T:0.01:T;

figure,
plot(t,X);xlabel('Time');ylabel('Amplitude');title('Input Function');

figure,
plot(abs(x2));
xlabel('frequency');ylabel('Amplitude');title('Magnitude Response of FRFT');

figure,
plot(imag(x2));
xlabel('frequency');ylabel('Amplitude');title('Imaginary part of FRFT');
 
toc