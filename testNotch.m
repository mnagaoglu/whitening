% test notch filtering 
clc
close all
clearvars;
Fs = 540;
t = 0:1/Fs:1;
fc = 60;
N = 100;
Ye = zeros(N,ceil(length(t)/2));
Yef = zeros(N,ceil(length(t)/2));
Yef2 = zeros(N,ceil(length(t)/2));
A = 100;
px2deg = 10/512;



% notch filter
d = designfilt('bandstopiir','FilterOrder',2, ...
       'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
       'DesignMethod','butter','SampleRate',Fs);

figure;cla;
for j=0:5:50
    for i=1:N
        % e = randn(length(t),1) + 0.1*cos(2*pi*t*fc)';
        e = round(A*(pinknoise(length(t)) + j*.01*cos(2*pi*(fc+1*randn(1,length(t))).*t)...
            + j*.01*cos(2*pi*(2*fc+1*randn(1,length(t))).*t))) * px2deg;

        ef = filtfilt(d,e);

        [~,~,Ye(i,:)] = getFFT(Fs,e);
        [~,~,Yef(i,:)] = getFFT(Fs,ef);
        [~,f,Yef2(i,:)] = getFFT(Fs,medfilt1(round(ef/px2deg),3)*px2deg);

    end


    loglog(f,mean(Ye.^2,1),'-k'); hold on;
    loglog(f,mean(Yef.^2,1),'-r'); hold on;
    loglog(f,mean(Yef2.^2,1),'-b'); hold on;
    xlabel('f (Hz)')
    ylabel('Power Spectrum')
    grid on;
end