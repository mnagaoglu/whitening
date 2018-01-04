function TestStitching

close all;
clc;

% iterations
N = 50;

% load the image, note that the original field of view for this image is 40
% deg in horizontal dimension.
fieldOfView = 40; % deg

% the following defines how long of a line will be taken from the image.
% This will effectively define the "width" of our space-time image.
lineLengthPx = 512; 

% define temporal parameters
T = 1; % seconds 
Fs = 512; % Hz
D = [100]; % arcmin^2/sec diffusion constant;
t = linspace(-T/2,T/2, T*Fs);

% get the list of images 
fList = dir([pwd filesep 'images' filesep '*.jpg']);

for di = 1:length(D)

    W = zeros(N,T*Fs);
    
    avgPSD = zeros(length(T*Fs), lineLengthPx);
    noMotionAvgPSD = zeros(lineLengthPx,1);
    
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

        % get a simulated eye movement trace in deg
        W(iter,:) = SimulateEyeMovements(T, Fs,D(di));

        % convert that to pixel units
        eyeMovementsPx = round(W(iter,:) / pixelSizeInDeg);        

        [x,y] = GetInitialPosition(height,width,eyeMovementsPx,lineLengthPx);

        % create the space-time image
        st = zeros(length(T*Fs), lineLengthPx);
        for i=1:length(t)
            st(i,:) = im(y, x+eyeMovementsPx(i) : x+eyeMovementsPx(i)+lineLengthPx-1);
        end
        
%         imwrite(st,'forest.jpg');
        
        % remove DC
        st = st - mean(st(:));

        % compute the frequency response
        [Fst,tf,sf] = GetFFT(t,lineLengthPx, st, pixelSizeInDeg);
        noFst = fftshift(abs(fft(im(y, x : x+lineLengthPx-1))).^2)';

        % sum up PSDs 
        avgPSD = avgPSD + abs(fftshift(Fst)).^2/(T*Fs*lineLengthPx);
        noMotionAvgPSD = noMotionAvgPSD + noFst/(T*Fs);
    end

    % normalize to make it average rather than sum
    avgPSD = avgPSD / N;
    noMotionAvgPSD = noMotionAvgPSD / N;

    figure('name',sprintf('D:%d',D(di)),'units',...
        'normalized','outerposition',[0.1891    0.0908    0.68    0.7375]); 
    subplot(2,3,1)
    plot(t,W','-r','LineWidth',2);
    hold on;
    plot(t,W(end,:)','-k','LineWidth',2);
    plot([0 0],[-.5 .5]*apertureSizeDeg + W(N,find(t>0,1))*[1 1],'-k','LineWidth',2);
    set(gca,'fontsize',16);
    xlabel('time (sec)');
    ylabel('position (deg)');
    ylim([min(W(N,:))-apertureSizeDeg/2 max(W(N,:))+apertureSizeDeg/2]);
    axis square;
    
    subplot(2,3,2)
    imagesc([min(W(N,:))-apertureSizeDeg/2 max(W(N,:))+apertureSizeDeg/2],[min(t) max(t)],st);
    colormap(gca,gray);
    ylabel('time (sec)');
    xlabel('position (deg)');
    set(gca,'fontsize',16);
    title('eye motion','color','r')
    axis square;
    
    subplot(2,3,3)
    imagesc([min(sf) max(sf)],[min(tf) max(tf)],(20*log10(avgPSD)));
    colormap(gca,gray);
    ylabel('temporal freq. (Hz)');
    xlabel('spatial freq. (cpd)');
    set(gca,'fontsize',16);
    axis square;
    

%     subplot(2,3,5)
%     plot(sf, 20*log10(sum(2*avgPSD, 1)),'-r','LineWidth',2);
%     hold on;
%     plot(sf, 20*log10(sum(2*noMotionAvgPSD  + 0.00001, 1)),'-b','LineWidth',2);
%     ylabel('power (dB)');
%     xlabel('spatial freq. (cpd)');
%     set(gca,'fontsize',16,'yscale','linear','xscale','log');
%     xlim([.3 60])
%     grid on;
    
    subplot(2,3,5)
    plot(sf, 20*log10(2*avgPSD(find(tf>10,1),:)),'-r','LineWidth',2);
    hold on;
    plot(sf, 20*log10(2*noMotionAvgPSD),'-b','LineWidth',2);
    ylabel('power (dB)');
    xlabel('spatial freq. (cpd)');
    set(gca,'fontsize',16,'yscale','linear','xscale','log');
    xlim([.3 60])
    grid on;
    axis square;
    
    subplot(2,3,6)
    plot(tf, 20*log10(sum(2*avgPSD, 2)),'-r','LineWidth',2);
    hold on;
    plot(tf, 20*log10(0.00001),'-b','LineWidth',2);
    ylabel('power (dB)');
    xlabel('temporal freq. (Hz)');
    set(gca,'fontsize',16,'yscale','linear','xscale','log');
    xlim([2 100])
    grid on;
    axis square;
end


function [x,y] = GetInitialPosition(height,width,eyeMovementsPx,lineLengthPx)

% get a random starting point on the image, avoiding the boundaries by
% taking into account the max and min eye position. we do not need to check
% boundaries for y, since we will be simulating only horizontal eye
% movements
y = randi(height,1); 
while true
    x = randi(width,1);
    if (x + min(eyeMovementsPx)) > 1 && (x + max(eyeMovementsPx) + lineLengthPx) < width
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
Fst = fft2(st);




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