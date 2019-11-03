% Additional_hammingwindow
% Zeyu Liu
% 10/13/2019

% draw hamming window
    %N = 640  
    %  w(n) = 0.54-0.46*cos(2*pi*n/N)
    %plot(w);

h = hamming(640);
plot(h);
pause;

%draw coincidence's hamming window 
% first Method
N=640;
t=1:N;
h1=hamming(N);
plot(t,h1);
hold on;    
plot(t+N/2,h1);
hold on;
plot(t+N,h1);
hold on;
plot(t+1.5*N,h1);
hold on;
plot(t+2*N,h1);
pause;

% Second Method
N=640;
t=1:N;
h1=hamming(N);

for i = 1:10
    plot(t+(i-1)*N/2,h1)
    hold on;
end;

