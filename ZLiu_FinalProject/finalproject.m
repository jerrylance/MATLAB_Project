% finalproject
% Zeyu Liu
% 11/29/2019
function finalproject
% 1.Read TIMIT data file into a vector x.
[x,Fs] = audioread('LDC93S1.wav');
N = length(x);
plot(x);
pause;
Fs = 16000;
% N = 46797, Fs = 16000
sound(x,Fs);
pause;

% 2.Apply u-Law compression to the data x that is quantized in 16 bits, 
% and convert to an 8-bit data y.

% uencode(x,16) Uniform quantization and encoding of the input into 16-bits.
% if not use the quantization, it will use more storage space, waste the
% storage source. In this file, uencode and udecode 16 step may omit.
a = uencode(x,16);
plot(a);
pause;
% udecode(x,16) Uniform decoding of the input.
b = udecode(a,16);
plot(b);
pause;
% compand() Source code mu-law or A-law compressor or expander.
c = compand(b,255,max(abs(b)),'mu/compressor');
plot(c);
pause;
% need quantization again, because mu/compress(ln() function, uneven) change the 
% interval size. But when we use ADPCM, we need the same interval of the y-axis.
a1 = uencode(c,8);
plot(a1);
pause;
y = udecode(a1,8);
plot(y);
pause;
% 3.Implement the IMA ADPCM algorithm.
% The function are adpcm_y = adpcm_encoder(raw_y) and raw_y = adpcm_decoder(adpcm_y)

% 4.Apply the ADPCM Encoder to compressed data y, and obtain encoded data z.
z = adpcm_encoder(y);
plot(z);
pause;
% 5.Apply the ADPCM Decoder to the encoded data z, and obtain decoded data y’.
y1 = adpcm_decoder(z);
plot(y1);
pause;

% 6.Apply u-Law expansion to the decoded data y’, 
% and obtain the reconstructed speech x’.
x1 = compand(y1,255,max(abs(y1)),'mu/expander');

% 7.Play and plot the reconstructed speech x’.
plot(x1);
pause;
sound(x1,Fs);
pause;
% wirte a file to the record result sound.
audiowrite('reconstructedx1.wav',x1,Fs);


% Additional work for analysis the u-Law compression
[x,Fs] = audioread('LDC93S1.wav');
N = length(x);
M = 320; % number of segments
Ns = floor(N/M); % number of samples per segments
for i = 1:50:M
    % start index = (i-1)*Ns+1, end index = i*Ns
    xs = x((i-1)*Ns+1:i*Ns);
    subplot(2,1,1),plot(xs);
    title('original data in different segments');
    pause;
    xs = uencode(xs,16);
    xs = udecode(xs,16);
    x2 = compand(xs,255,max(abs(b)),'mu/compressor');
    subplot(2,1,2),plot(x2);
    title('u-Law compression in different segments');
    pause;
end

% standard u-Law compression principle
xu = linspace(-1,1);
yu = xu;
u1 = compand(yu,255,max(abs(yu)),'mu/compressor');
u2 = compand(yu,63,max(abs(yu)),'mu/compressor');
u3 = compand(yu,7,max(abs(yu)),'mu/compressor');

subplot(1,1,1),plot(xu,u1,xu,u2,xu,u3,xu,yu);
legend('u = 255','u = 63','u = 7','y = [-1,1]','location','nw');
title('standard u-Law compression with different u');
grid on
pause;
% u-Law compression for the range of y1
xu = linspace(-0.1,0.1);
yu = xu;
u1 = compand(yu,255,max(abs(y1)),'mu/compressor');
u2 = compand(yu,63,max(abs(y1)),'mu/compressor');
u3 = compand(yu,7,max(abs(y1)),'mu/compressor');

subplot(1,1,1),plot(xu,u1,xu,u2,xu,u3,xu,yu);
legend('u = 255','u = 63','u = 7','y = [-0.1,0.1]','location','nw');
title('u-Law compression with different u for range of y1');
grid on
pause;

% A example for quantization noise
u = linspace(-1,1);
uu = uencode(u,3);
uuu = udecode(uu,3);
uu1 = uencode(u,5);
uuu1 = udecode(uu1,5);
plot(u,uuu);
hold on
plot(u,uuu1);
hold on
plot(u,u);
title('A example for quantization noise');
legend('quantized in 3 bits','quantized in 5 bits','original signal','location','nw');