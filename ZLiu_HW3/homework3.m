% Homework 3
% Window operations
% Zeyu Liu
% 10/8/2019

function homework3
% 1(a). window of length L.(L = 401)
L = 401;
h = zeros(L);
for n=1:L,
    h(n) = 0.54 - 0.46*cos(2*pi*n/L);
end;
plot(h);
pause;

%or use tese is the same
%h1 = hamming(L);
%plot(h1);

% 1(b). Apply FFT
% h = hamming(L);
% plot(h1);
% pause;
h1 = double(hamming(L));
plot(h1);
pause;

% H = fft(h);
% logh = log(abs(H))+1;
% plot(logh);
% pause;
% plot(angle(H));
% pause;

H1 = fft(h1);
% plot magnitude
logh = log(abs(H1))+1;
plot(logh);
pause;
% plot phase
plot(angle(H1));
pause;

% 2. compute short-time energy and apply to the TIMIT speech sample
fname = 'LDC93s1.wav';
x = audioread(fname);
% short-time energy E = convolution(x^2(n)*h(n)), h(n) = w^2(n)
L = [51, 101, 201, 401];
for k = 1:4,
    w = hamming(L(k));
    h = w.^2;
    y = conv(x.^2,h);
    N = length(y);
    for m = 1:floor(N/4),
        energy(m) = y((m-1)*4+1);
    end;
    plot(x);
    pause;
    plot(energy);
    pause;
end;

% 3. short-time zero-crossings
for k = 1:4,
    w = hamming(L(k));
    h = w.^2;
    N = length(x)-1;
    diff = zeros(1,N);
    for n = 1:N,
        diff(n) = abs(sign(x(n+1))-sign(x(n)));
    end;
    y = conv(diff,h);
    M = length(y);
    for m = 1:floor(M/4);
        zc(m) = y((m-1)*4+1);
    end;
    plot(x);
    pause;
    plot(zc);
    pause;
end;

% 4.compute short-time auto-correlation function
% test with TIMIT speech sample.  Display results 
% as a 2D image (time lag k as y axis, short-time 
% index as x axis).

% short-time auto-correlation function:
% R = convolution(x(n)*x(n-k)) or convolution(x(n+k)*x(n))
fname = 'LDC93s1.wav';
x = audioread(fname);
L = [51, 101, 201, 401];
% time lag k as y axis, short-time index as x axis
for k1 = 1:4,
    w = hamming(L(k1));
    N = length(x)-640;
    co = zeros(floor(N/320),250)
    for i=1:320:N,
        for k=1:250,
            x1 = x(i:i+320);
            y1 = conv(x1,w);
            x2 = x(i+k:i+k+320);
            y2 = conv(x2,w);
            tk(k) = y1'*y2;
        end;
        co((i-1)/320+1,:) = tk; % 
    end;    
    imagesc(co');
    colormap(gray);
    pause;
end;
