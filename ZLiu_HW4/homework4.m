% Homework 4
% Zeyu Liu
% 11/12/2019
function homework4
s = rng(1);% use a seed control the result
p = 20;
N = 1024;

% 1. generate white Gaussian random noise w
w = randn(N,1);
plot(w);
pause;
% 2. generate the 1024-point signal x 
x = zeros(N,1);
for n = 1:N,
    x(n) = cos(0.1*pi*n) + 0.2*sin(0.2*pi*n) + 0.2*w(n)
end;
plot(x);
pause;
% 3. Estimate the parameters of 20th order autoregressive (linear prediction) model by a direct inversion method.
zz = zeros(p);
zy = zeros(p,1);
for n = p+1:N,
    z = x(n-p:n-1);
    zz = zz + z*z';
    zy = zy + z*x(n);
end;

% estimate by inversion
a = inv(zz)*zy;
a = flip(a); % returns with the order of the elements reversed. 
a = [1; -a];
plot(a);
pause;
% 4. Compute the power spectral density using the estimated parameters, and plot the spectrum
phi = zeros(p+1,1);
s = zeros(N,1);
for n = 1:N,
    omega = (n-1)*pi/N;
    for m = 1:p+1,
        phi(m) = exp(1j*(m-1)*omega);
    end;
    S(n) = 1/abs(a'*phi)^2;
end;
plot(log(S));
pause;
% The peak in the plot is x1 = 103, and x2 = 207
omega1 = (103-1)*pi/N % omega1 = 0.3129
omega2 = (207-1)*pi/N % omega2 = 0.6320
pause;
% 5. Find the roots of A(z)
r = roots(a);
r
% can also use the pzmap(),the image is same
% pzmap(1,r);
plot(r,'x','LineWidth',3,'MarkerSize',8);
hold on;
viscircles([0,0],1,'LineWidth',1);
hold off;
pause;

% find the spectral peaks from the computed roots.
% The largest two peaks in part 5 are the point nearest from the unit circle.
% Caculate the theta of the two points, and compare with the omega in part 4.

Re = real(r);
Im = imag(r);
theta = atan(Im./Re);
theta
% theta1 = 0.3143
% theta2 = 0.6328
pause;

% 6. Yule-Walker
a = aryule(x,p);
plot(a);
pause;
% Compute the power spectral density
phi = zeros(p+1,1);
s = zeros(N,1);
for n = 1:N,
    omega = (n-1)*pi/N;
    for m = 1:p+1,
        phi(m) = exp(1j*(m-1)*omega);
    end;
    S(n) = 1/abs(a*phi)^2;
end;
plot(log(S));
pause;

% The peak in the plot is x1 = 103, and x2 = 207
omega1 = (103-1)*pi/N % omega1 = 0.3129
omega2 = (207-1)*pi/N % omega2 = 0.6320
pause;

% Find the roots 
r = roots(a);
r
plot(r,'x','LineWidth',3,'MarkerSize',8);
hold on;
viscircles([0,0],1,'LineWidth',1);
hold off;
pause;

% find the spectral peaks from the computed roots.
Re = real(r);
Im = imag(r);
theta = atan(Im./Re);
theta
% theta1 = 0.3141
% theta2 = 0.6325
pause;

% 7. Burg-Anderson 
a = arburg(x,p);
plot(a);
pause;
% Compute the power spectral density
phi = zeros(p+1,1);
s = zeros(N,1);
for n = 1:N,
    omega = (n-1)*pi/N;
    for m = 1:p+1,
        phi(m) = exp(1j*(m-1)*omega);
    end;
    S(n) = 1/abs(a*phi)^2;
end;
plot(log(S));
pause;

% The peak in the plot is x1 = 103, and x2 = 207
omega1 = (103-1)*pi/N % omega1 = 0.3129
omega2 = (207-1)*pi/N % omega2 = 0.6320
pause;

% Find the roots 
r = roots(a);
r
plot(r,'x','LineWidth',3,'MarkerSize',8);
hold on;
viscircles([0,0],1,'LineWidth',1);
hold off;
pause;

% find the spectral peaks from the computed roots.
Re = real(r);
Im = imag(r);
theta = atan(Im./Re);
theta
% theta1 = 0.3143
% theta2 = 0.6328
pause;

