function TestIn1D
% This function is written to perform the spectral analyses in 1D space and
% time. It is a better playground for medium to low-end computers and 
% computers without discrete GPUs. 
%
% To get realistic and naturalistic content, 100 images from 
% http://natural-scenes.cps.utexas.edu/db.shtml were used (~40deg field of 
% view horizontally).
%
%


% close all force;
clc;

% to use real eye movement traces instead of simulating them as Brownian
% motion, set the following flag to 1.
isRealEM = 1;
if isRealEM
    global organizedData;
    filePath = '../../ganglion/organizedData.mat';
    load(filePath,'organizedData');
    groupNo = 2; % 1:age-matched, 2:young-control, 3:AMD
    isFilter = 0; % eye movement traces were already filtered
    
    if isFilter
        % lowpass filtering eye positions 
        global lpFilt;
        global notchFilt;
        lpFilt = designfilt('lowpassfir','SampleRate',480,...
             'PassbandFrequency',120, ...
             'StopbandFrequency',150,'PassbandRipple',0.05, ...
             'StopbandAttenuation',60,'DesignMethod','kaiserwin');

        % notch filter
        notchFilt = designfilt('bandstopiir','FilterOrder',2, ...
                   'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
                   'DesignMethod','butter','SampleRate',480);
    end
end

is2D = 1; % set to 0, if you want to fix vertical position and just use horizontal

% iterations
N = 100;

% load the image, note that the original field of view for this image is 40
% deg in horizontal dimension.
fieldOfView = 40; % deg

% the following defines how long of a line will be taken from the image.
% This will effectively define the "width" of our space-time image.
lineLengthPx = 512; 

% define temporal parameters
T = 1; % seconds 
Fs = 512; % Hz
D = 100; % arcmin^2/sec diffusion constant;
% Dratios = [.01 .1 .2 .5 1 2 5 10 100];
Dratios = 1;
% actualDx = sqrt(D^2./(Dratios.^2 + 1));
actualDx = (D.*Dratios./(Dratios + 1));
actualDy = actualDx ./ Dratios;
t = linspace(-T/2,T/2, T*Fs);

forLegend = [];
for i=1:length(Dratios)
    forLegend{i} = mat2str(Dratios(i));
end

% get the list of images 
fList = dir([pwd filesep 'images' filesep '*.jpg']);


wb = waitbar(0,'Please wait...');

for di = 1:length(actualDx)

    Wx = zeros(N,T*Fs);
    Wy = zeros(N,T*Fs);
    
    avgPSD = zeros(length(T*Fs), lineLengthPx);
    noMotionAvgPSD = avgPSD;
    
    
    for iter = 1:N
        
        % get a random image
        imageID = randi(length(fList),1);
        im = imread([fList(imageID).folder filesep fList(imageID).name]);
        im = double(im)/255;
        
        % get the size of the image in pixels
        [height, width] = size(im);

        % compute the pixel size in deg
        pixelSizeInDeg = fieldOfView / width;

        % find the corresponding "aperture size" for the selected line width
        apertureSizeDeg = pixelSizeInDeg * lineLengthPx;

        % get EM traces
        if isRealEM
            [Wx(iter,:), Wy(iter,:), t] = GetRealEyeMovements(T,Fs,groupNo, isFilter);
        else
            % get a simulated eye movement trace in deg
            Wx(iter,:) = SimulateEyeMovements(T, Fs,actualDx(di));
            Wy(iter,:) = SimulateEyeMovements(T, Fs,actualDy(di));
        end

        % convert that to pixel units
        eyeMovementsPxX = round(Wx(iter,:) / pixelSizeInDeg);     
        eyeMovementsPxY = round(Wy(iter,:) / pixelSizeInDeg);

        [x,y] = GetInitialPosition(height,width,eyeMovementsPxX,eyeMovementsPxY,lineLengthPx);

        % create the space-time image
        st = zeros(length(T*Fs), lineLengthPx);
        for i=1:length(t)
            if is2D
                st(i,:) = im(y+eyeMovementsPxY(i) , ...
                             x+eyeMovementsPxX(i) : x+eyeMovementsPxX(i)+lineLengthPx-1);
            else
                st(i,:) = im(y, ...
                             x+eyeMovementsPxY(i) : x+eyeMovementsPxY(i)+lineLengthPx-1);
            end
        end

        % compute the frequency response
        [Fst,tf,sf] = GetFFT(t,lineLengthPx, st - mean(st(:)), pixelSizeInDeg);
        noFst = abs(fftshift(fft(st - repmat(mean(st,2),1,size(st,2)) ,[],2)));

        % sum up PSDs 
        avgPSD = avgPSD + (2*Fst).^2/(T*Fs*lineLengthPx);
        noMotionAvgPSD = noMotionAvgPSD + (2*noFst).^2/(T*Fs);
        
        
        waitbar(iter/N,wb,sprintf('Dx: %d, Dy: %d, progress %% %0.f', ...
            round(actualDx(di)), round(actualDy(di)), 100*iter/N))
        
    end

    % normalize to make it average rather than sum
    avgPSD = avgPSD / N;
    noMotionAvgPSD = noMotionAvgPSD / N;

    if isRealEM
        titleText = 'With real eye movements';
    else
        titleText = sprintf('Dx: %d, Dy: %d', ...
            round(actualDx(di)), round(actualDy(di)));
    end
    
    if di == 1
        fh = figure('name',titleText,'units',...
            'normalized','outerposition',[-1.2607    0.5781    1.2607    0.4000]); 
        subplot(1,5,1)
        plot(t,Wx','-r','LineWidth',2);
        hold on;
        plot(t,Wx(end,:)','-k','LineWidth',2);
    %     plot([0 0],[-.5 .5]*apertureSizeDeg + Wx(N,find(t>0,1))*[1 1],'-k','LineWidth',2);
        set(gca,'fontsize',16);
        xlabel('time (sec)');
        ylabel('position (deg)');
    %     ylim([min(Wx(N,:))-apertureSizeDeg/2 max(Wx(N,:))+apertureSizeDeg/2]);
        axis square;
        box off

        subplot(1,5,2)
        imagesc([min(Wx(N,:))-apertureSizeDeg/2 max(Wx(N,:))+apertureSizeDeg/2],[min(t) max(t)],st);
        colormap(gca,gray);
        ylabel('time (sec)');
        xlabel('position (deg)');
        set(gca,'fontsize',16);
        axis square;
        box off;

        subplot(1,5,3)
        imagesc([min(sf) max(sf)],[min(tf) max(tf)],(10*log10(avgPSD)));
        colormap(gca,gray);
        ylabel('temporal freq. (Hz)');
        xlabel('spatial freq. (cpd)');
        set(gca,'fontsize',16);
        axis square;
        box off;
    end
    
    figure(fh);
    subplot(1,5,4)
    if di == 1
        plot(sf, mean(10*log10(noMotionAvgPSD),1),'-b','LineWidth',2);
    end
    hold on;
    
%     cols = jet(T*Fs/2-1);
%     for lx = 1 : T*Fs/2-1
%         plot(sf, 10*log10(avgPSD(T*Fs/2 + lx,:)),'-','Color',cols(lx,:),'LineWidth',2);
%         hold on;
%         plot(sf, 10*log10(noMotionAvgPSD(T*Fs/2 + lx,:)),'--','Color',cols(lx,:),'LineWidth',2);
%     end
    plot(sf, mean(10*log10(avgPSD),1),'-','Color',[1 0 0]*di/length(actualDx),'LineWidth',2);
    hold on;
    ylabel('power (dB)');
    xlabel('spatial freq. (cpd)');
    set(gca,'fontsize',16,'yscale','linear','xscale','log');
    xlim([.5 60])
    grid on;
    axis square;
    box off;

    
    subplot(1,5,5)
    plot(tf, mean(10*log10(avgPSD), 2),'-','Color',[1 0 0]*di/length(actualDx),'LineWidth',2);
    hold on;
%     plot(tf, 10*log10(0.00001),'-b','LineWidth',2);
    ylabel('power (dB)');
    xlabel('temporal freq. (Hz)');
    set(gca,'fontsize',16,'yscale','linear','xscale','log');
    xlim([2 100])
    grid on;
    axis square;
    box off;
    if di == length(actualDx)
        lh = legend(forLegend,'location','best');
        lh.Title.String = 'Diffusion ratio x/y';
    end
    
    % break if real EM was used since D is irrelevant in this context
    if isRealEM
        break;
    end
    
    
end

delete(wb);


function [x,y] = GetInitialPosition(height,width,eyeMovementsPxX,eyeMovementsPxY,lineLengthPx)

% get a random starting point on the image, avoiding the boundaries by
% taking into account the max and min eye position. 
while true
    x = randi(width,1);
    if (x + min(eyeMovementsPxX)) > 1 && (x + max(eyeMovementsPxX) + lineLengthPx) < width
        break;
    end
end


while true
    y = randi(height,1); 
    if (y + min(eyeMovementsPxY)) > 1 && (y + max(eyeMovementsPxY) + lineLengthPx) < height
        break;
    end
end




function [Fst,tf,sf] = GetFFT(t,lineLengthPx, st, pixelSizeInDeg)

% compute frequency axes
Lt = length(t);
Ls = lineLengthPx;
dt = diff(t(1:2));
tf = (1/dt) * linspace(-.5, .5, Lt);
sf = (1/pixelSizeInDeg) * linspace(-.5, .5, Ls);

% take the FFt
Fst = abs(fftshift(fft2(st)));




function [W, t, dt] = SimulateEyeMovements(T,Fs,D)

% Definition of variables
if nargin<2
    T = 1;  % second
    Fs = 512; % number of samples per second
end

diffusion = D/3600; % deg^2/sec
nsamples = 1; % indicates the number of instances of Brownian motion 
t = linspace(-T/2,T/2, T*Fs); % time array for plotting
dt = diff(t(1:2));

% Standard Brownian motion (with correlation)
dW = sqrt(2*diffusion*dt) * randn(nsamples,length(t)-1); 
W = [zeros(nsamples,1) cumsum(dW,2)];



function [Wx, Wy, t, dt] = GetRealEyeMovements(T,Fs, groupNo, isFilter)

global organizedData;

% normal young observers, groupNo=2
groups = [organizedData.group];
youngIndices = find(groups == groupNo);

% get a random subject
subjectIndex = randi(length(youngIndices),1);

% get the time and position traces
time = organizedData(youngIndices(subjectIndex)).stitchedTime;
pos = organizedData(youngIndices(subjectIndex)).stitchedPosition;

if isFilter
    pos = FilterEM(pos);
end

% how many samples requested
nSamples = T*Fs;

% get a random chunk of data
actualDt = diff(time(1:2));
actualSamples = round(T/actualDt);
startIndex = randi(length(time)-actualSamples-1,1);

% interpolate time and pos according to the requested sampling rate

newTime = (0:1/Fs:(1/Fs)*(nSamples-1))' + time(startIndex);
newPos(:,1) = interp1(time(startIndex:startIndex+actualSamples-1), ...
                pos(startIndex:startIndex+actualSamples-1,1), ...
                newTime, 'pchip');
newPos(:,2) = interp1(time(startIndex:startIndex+actualSamples-1), ...
                pos(startIndex:startIndex+actualSamples-1,2), ...
                newTime, 'pchip');
            
dt = 1/Fs;
W = newPos - repmat(newPos(1,:),length(newTime),1);
t = newTime - newTime(1);

Wx = W(:,1);
Wy = W(:,2);

% figure(1000);
% cla;
% plot(time(startIndex:startIndex+actualSamples-1)-time(startIndex),...
%     pos(startIndex:startIndex+actualSamples-1,:) - pos(startIndex,:),...
%     '-r','LineWidth',2); 
% hold on;
% plot(t,W,'-b','LineWidth',2);
% xlabel('time (sec)')



function newPos = FilterEM(pos)

global lpFilt;
global notchFilt;

% median filt
w = 7;
mfposition(:,1) = medfilt1(pos(:,1),w);
mfposition(:,2) = medfilt1(pos(:,2),w);

% lowpass filt
lpmfposition(:,1) = filtfilt(lpFilt,mfposition(:,1));
lpmfposition(:,2) = filtfilt(lpFilt,mfposition(:,2));

% notch filt
newPos(:,1) = filtfilt(notchFilt,lpmfposition(:,1));
newPos(:,2) = filtfilt(notchFilt,lpmfposition(:,2));

% figure(100);
% cla;
% plot(pos); hold on;
% plot(newPos,'LineWidth',2);
% title('filtered signal')


