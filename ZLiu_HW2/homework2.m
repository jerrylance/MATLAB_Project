% Homework 2
% Exercise in reading and playing speech files
% and computing spectrogram
% Zeyu Liu
% 10/8/2019

function homework2

% Write a matlab code to read the wave file into a vector x.
[x,Fs] = audioread('LDC93S1.wav');
N = length(x);
plot(x);
pause;
% N
Fs = 16000;
% N = 46797, Fs = 16000, N/5 = 9359.4
sound(x,Fs);
pause;

% Divide the vector x into five segments and plot them.
M = 5; % number of segments
Ns = floor(N/M); % number of samples per segments
for i = 1:M,
    % start index = (i-1)*Ns+1, end index = i*Ns
    xs = x((i-1)*Ns+1:i*Ns);
    plot(xs);
    sound(xs,Fs);
    pause;
end;

% Divide the speech file into 20ms segments (320 samples each)
% compute periodogram

Ns = 320 % number of samples per segments
M = floor(N/320); % number of segments
B = zeros(160,M); % blank image, 160 = 320/2, because symmetrical
% M = 146
for i = 1:M,
    % start index = (i-1)*Ns+1, end index = i*Ns
    xs = x((i-1)*Ns+1:i*Ns);
    plot(xs);
    %sound(xs,Fs);
    pause;
    % compute spectogram
    Xs = fft(xs);
    fx = (abs(Xs)).^2/Ns;
    plot(fx);
    pause;
    plot(log(fx));
    pause;
    % Then compute spectrogram by constructing an image made of set of 
    %periodograms. 
    % Display the negative image in log-scale. 
    B(:,i) = log(fx(1:160));
end;

imagesc(B);
colormap(gray);
pause;
imagesc(-B);
colormap(gray);
pause;




