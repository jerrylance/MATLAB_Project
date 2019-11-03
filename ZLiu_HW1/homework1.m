% Homework 1
% Generate a signal of two sinusoids with additive white noise
% then plot, apply DFT, and computer power spectrum.
% Zeyu Liu
% Created 9/17/2019
% Revised

function homework1
% 1. Generate 1024-point white Gaussian random noise w, 
% with zero mean and unit variance.

w = randn(1024,1);
plot(w);
pause;

% 2. Generate the 1024-point signal x of two sinusoids with noise 
x = zeros(1024,1);

for n= 1:1024;
    x(n) = cos(0.1*pi*n) + 0.2*sin(0.2*pi*n) + 0.2*w(n);
end;

% 3. plot x(n)
plot(x);
pause;

% 4. Compute the DFT of x(n), and obtain X(k).
X = fft(x);
plot(X);
pause;
% plot real plot
plot(real(X));
pause;
% plot imaginary part
plot(imag(X));
pause;
% plot magnitude
plot(abs(X));
pause;
% plot phase
plot(angle(X));
pause;

% 5. Plot the periodogram, for k=1,..,512. How many peaks can you observe? What are the locations of the peaks?

for m = 1:512
    X1(m) = (abs(X(m)))^2;
end;
plot(X1);
pause;
% plot real plot
% plot(real(X1));
% pause;
% plot imaginary part
% plot(imag(X1));
% pause;
% plot magnitude
% plot(abs(X1));
% pause;
% plot phase
% plot(angle(X1));
% pause;

% 6. Repeat the above experiment by changing the frequency of the second sinusoid. Do you observe any differences?

for n= 1:1024,
    x2(n) = cos(0.1*pi*n) + 0.2*sin(0.12*pi*n) + 0.2*w(n);
end;
plot(x2);
pause;
X2 = fft(x2);
plot(X2);
pause;
% plot real plot
plot(real(X2));
pause;
% plot imaginary part
plot(imag(X2));
pause;
% plot magnitude
plot(abs(X2));
pause;
% plot phase
plot(angle(X2));
pause;

for m1 = 1:512
    X3(m1) = (abs(X2(m1)))^2;
end;
plot(X3);
pause;
