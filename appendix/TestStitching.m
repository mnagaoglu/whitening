function TestStitching
% test to see if using a longer eye motion trace in time will change power
% spectra...
% hint, it does not, except for the change in overal energy, which is
% expected.

close all;
clc;

isPlotTime = 0;

% iterations
N = 5000;

% define temporal parameters
T = fliplr([1/8 1/4 1/2 1 2]); % seconds 
Fs = 512; % Hz
D = 40; % arcmin^2/sec diffusion constant;

figure('units','normalized','outerposition',[-0.6750   -0.2000    0.6750    0.5000]);
for i=1:length(T)

    FW = zeros(1,round(T(i)*Fs));

    try
        for iter=1:N
            % simulate eye motion
            [W, t, dt] = SimulateEyeMovements(T(i),Fs,D);

            % take fft
            tempFW = abs(fftshift(fft(W)));
            FW = FW + 2* (tempFW.^2 / (T(i)*Fs));

            if isPlotTime
                subplot(1,2,1)
                plot(t, W,'-','LineWidth',2,'Color',[1 0 0]*i/length(T));
                hold on;
                ylabel('position (deg)');
                xlabel('time (sec)');
                set(gca,'fontsize',16,'yscale','linear','xscale','linear');
                grid on;
                axis square;

            end
        end
    catch err
        rethrow(err);
    end

    % take average
    tf = (1/dt) * linspace(-.5, .5, T(i)*Fs);
    FW = FW / N;

    
    
    subplot(1,2,2)
    plot(tf, 20*log10(FW),'-','LineWidth',2,'Color',[1 0 0]*i/length(T));
    hold on;
    ylabel('power (dB)');
    xlabel('temporal freq. (Hz)');
    set(gca,'fontsize',16,'yscale','linear','xscale','log');
    xlim([1 100])
    grid on;
    axis square;
    
    
end








function [W, t, dt] = SimulateEyeMovements(T,Fs,D)

% Definition of variables
if nargin<2
    T = 1;  % second
    Fs = 512; % number of samples per second
end

diffusion = D/3600; % deg^2/sec
nsamples = 1; % indicates the number of instances of Brownian motion 
t = linspace(0,T, T*Fs); % time array for plotting
dt = diff(t(1:2));

% Standard Brownian motion (with correlation)
dW = sqrt(2*diffusion*dt) * randn(nsamples,length(t)-1); 
W = [zeros(nsamples,1) cumsum(dW,2)];




