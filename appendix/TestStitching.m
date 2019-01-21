function TestStitching
% test to see if using a longer eye motion trace in time will change power
% spectra...
% hint, it does not, except for the change in overal energy, which is
% expected.
%
%

close all;
clc;

isRealEM = 0;


if isRealEM
    global organizedData;
    filePath = '../../ganglion/organizedData.mat';
    load(filePath,'organizedData');
end

isPlotTime = 1;
simulateFs = 1;

% iterations
N = 100;

% define temporal parameters
if ~simulateFs
    T = fliplr([1/16 1/8 1/4 1/2 1 2 4]); % seconds 
    Fs = repmat(1024,1,length(T)); % Hz
else
    Fs = fliplr([128 256 512 1024 2048]); % Hz
    T = ones(1,length(Fs)); % seconds 
end

D = 40; % arcmin^2/sec diffusion constant;

figure('units','normalized','outerposition',[-0.6750   -0.2000    0.5554    0.3895]);
for i=1:length(T)

    FW = zeros(1,round(T(i)*Fs(i)));

    try
        for iter=1:N
            
            if ~isRealEM
                % simulate eye motion
                [W, t, dt] = SimulateEyeMovements(T(i),Fs(i),D);
            else
                % get real eye motion
                [W, t, dt] = GetRealEyeMovements(T(i),Fs(i));
            end
            
            % take fft
            tempFW = 2*abs(fftshift(fft(W)));
            FW = FW + (tempFW.^2 / (T(i)*Fs(i)));

            if isPlotTime
                subplot(1,2,1)
                plot(t, W,'-','LineWidth',1,'Color',[1 0 0]*i/length(T));
                hold on;
                ylabel('position (deg)');
                xlabel('time (sec)');
                set(gca,'fontsize',20,'yscale','linear','xscale','linear');
                grid on;
                axis square;
                box off;

            end
        end
    catch err
        rethrow(err);
    end

    % take average
    tf = (1/dt) * linspace(-.5, .5, T(i)*Fs(i));
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

% put a legend
leg = [];
for i=1:length(T)
    if simulateFs
        leg{i} = mat2str(Fs(i));
        if i==1
            title('Sampling rate');
        end
    else
        if i==1
            title('Time window');
        end
        leg{i} = mat2str(T(i));
    end
end
legend(leg);




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

W = W - W(1);




function [W, t, dt] = GetRealEyeMovements(T,Fs)

global organizedData;

% normal young observers, group=2
groups = [organizedData.group];
youngIndices = find(groups == 2);

% get a random subject
subjectIndex = randi(length(youngIndices),1);

% get the time and position traces
time = organizedData(youngIndices(subjectIndex)).stitchedTime;
pos = organizedData(youngIndices(subjectIndex)).stitchedPosition;

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








