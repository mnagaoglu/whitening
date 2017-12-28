function [Y,f,P1] = getFFT(Fs,position)

L = length(position); % Length of signal

Y = fft(position);
P2 = abs(Y/L);
P1 = P2(1:floor(L/2+1));
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;

