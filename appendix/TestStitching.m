function TestStitching
% test to see if using a longer eye motion trace in time will change power
% spectra...
% hint, it does not, except for the change in overal energy, which is
% expected.
%
%

close all;
clc;

global organizedData;
filePath = '../../ganglion/organizedData.mat';
load(filePath,'organizedData');


isPlotTime = 0;

% iterations
N = 5000;

% define temporal parameters
T = fliplr([1/8 1/4 1/2 1 2 4]); % seconds 
Fs = 512; % Hz
D = 40; % arcmin^2/sec diffusion constant;

figure('units','normalized','outerposition',[-0.6750   -0.2000    0.6750    0.5000]);
for i=1:length(T)

    FW = zeros(1,round(T(i)*Fs));

    try
        for iter=1:N
            
%             % simulate eye motion
%             [W, t, dt] = SimulateEyeMovements(T(i),Fs,D);
            
            % get real eye motion
            [W, t, dt] = GetRealEyeMovements(T(i),Fs);

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
t = linspace(0,T, T*Fs); % time array for plotting
dt = diff(t(1:2));

% Standard Brownian motion (with correlation)
dW = sqrt(2*diffusion*dt) * randn(1,length(t)-1); 
W = [zeros(1,1) cumsum(dW,2)];




function [W, t, dt] = GetRealEyeMovements(T,Fs)

global organizedData;

% normal young observers, group=2
groups = [organizedData.group];
youngIndices = find(groups == 2);

% get a random subject
subjectIndex = randi(length(youngIndices),1);

% get the time and position traces
time = organizedData(youngIndices(subjectIndex)).newTime;
pos = organizedData(youngIndices(subjectIndex)).newPos;

% how many samples requested
nSamples = T*Fs;

% get a random chunk of data
actualDt = diff(time(1:2));
actualSamples = round(T/actualDt);
startIndex = randi(length(time)-actualSamples-1,1);

% interpolate time and pos according to the requested sampling rate

newTime = interp1(0:actualDt:actualDt*(actualSamples-1), ...
                time(startIndex:startIndex+actualSamples-1), ...
                0:1/Fs:(1/Fs)*(nSamples-1), 'pchip');
newPos = interp1(time(startIndex:startIndex+actualSamples-1), ...
                pos(startIndex:startIndex+actualSamples-1,2), ...
                newTime, 'pchip');
            
dt = 1/Fs;
W = newPos - newPos(1);
t = newTime - newTime(1);

% figure;
% plot(time(startIndex:startIndex+actualSamples-1)-time(startIndex),...
%     pos(startIndex:startIndex+actualSamples-1,2) - pos(startIndex,2),...
%     '-r','LineWidth',2); 
% hold on;
% plot(t,W,'-b','LineWidth',2);
% xlabel('time (sec)')








