function TestStitching

close all;
clc;

N = 1000;

% load the image, note that the original field of view for this image is 40
% deg in horizontal dimension.
fieldOfView = 40; % deg
im = imread('cps201004281302.jpg');
im = double(im)/255;

% get the size of the image in pixels
[height, width] = size(im);

% compute the pixel size in deg
pixelSizeInDeg = fieldOfView / width;

% the following defines how long of a line will be taken from the image.
% This will effectively define the "width" of our space-time image.
lineLengthPx = 512; 

% find the corresponding "aperture size" for the selected line width
apertureSizeDeg = pixelSizeInDeg * lineLengthPx;

% define temporal parameters
T = 1; % seconds 
Fs = 256; % Hz
D = [100]; % arcmin^2/sec diffusion constant;
t = linspace(-T/2,T/2, T*Fs);


for di = 1:length(D)

    W = zeros(N,T*Fs);
    
    avgPSD = zeros(length(T*Fs), lineLengthPx);
    noMotionAvgPSD = avgPSD;
    
    for iter = 1:N

        % get a simulated eye movement trace in deg
        W(iter,:) = SimulateEyeMovements(T, Fs,D(di));

        % convert that to pixel units
        eyeMovementsPx = round(W(iter,:) / pixelSizeInDeg);        

        [x,y] = GetInitialPosition(height,width,eyeMovementsPx,lineLengthPx);

        % create the space-time image
        st = zeros(length(T*Fs), lineLengthPx);
        noMotionSt = st;
        for i=1:length(t)
            st(i,:) = im(y, x+eyeMovementsPx(i) : x+eyeMovementsPx(i)+lineLengthPx-1);
            noMotionSt(i,:) = im(y, x : x+lineLengthPx-1);
        end
        
        % remove DC
        st = st - mean(st(:));
        noMotionSt = noMotionSt - mean(noMotionSt(:));

        % compute the frequency response
        [Fst,tf,sf] = GetFFT(t,lineLengthPx, st, pixelSizeInDeg);
        noFst = GetFFT(t,lineLengthPx, noMotionSt, pixelSizeInDeg);

        % sum up PSDs 
        avgPSD = avgPSD + abs(fftshift(Fst)).^2/(T*Fs*lineLengthPx);
        noMotionAvgPSD = noMotionAvgPSD + abs(fftshift(noFst)).^2/(T*Fs*lineLengthPx);
    end

    % normalize to make it average rather than sum
    avgPSD = avgPSD / N;
    noMotionAvgPSD = noMotionAvgPSD / N;

    figure('name',sprintf('D:%d',D(di)),'units',...
        'normalized','outerposition',[-2 -1 2 0.4]); 
    subplot(1,7,1)
    plot(t,W','-r','LineWidth',2);
    hold on;
    plot(t,W(end,:)','-k','LineWidth',2);
    plot([0 0],[-.5 .5]*apertureSizeDeg + W(N,find(t>0,1))*[1 1],'-k','LineWidth',2);
    set(gca,'fontsize',16);
    xlabel('time (sec)');
    ylabel('position (deg)');
    ylim([min(W(N,:))-apertureSizeDeg/2 max(W(N,:))+apertureSizeDeg/2]);

    subplot(1,7,2)
    imagesc([min(W(N,:))-apertureSizeDeg/2 max(W(N,:))+apertureSizeDeg/2],[min(t) max(t)],st);
    colormap(gca,gray);
    ylabel('time (sec)');
    xlabel('position (deg)');
    set(gca,'fontsize',16);
    title('eye motion','color','r')

    subplot(1,7,3)
    imagesc([min(W(N,:))-apertureSizeDeg/2 max(W(N,:))+apertureSizeDeg/2],[min(t) max(t)],noMotionSt);
    colormap(gca,gray);
    ylabel('time (sec)');
    xlabel('position (deg)');
    set(gca,'fontsize',16);    
    title('no motion','color','b')
    
    subplot(1,7,4)
    imagesc([min(sf) max(sf)],[min(tf) max(tf)],(log10(avgPSD)));
    colormap(gca,gray);
    ylabel('temporal freq. (Hz)');
    xlabel('spatial freq. (cpd)');
    set(gca,'fontsize',16);
    
    
    subplot(1,7,5)
    imagesc([min(sf) max(sf)],[min(tf) max(tf)],(log10(noMotionAvgPSD)));
    colormap(gca,gray);
    ylabel('temporal freq. (Hz)');
    xlabel('spatial freq. (cpd)');
    set(gca,'fontsize',16);

    subplot(1,7,6)
    plot(sf, mean((20*log10(2*avgPSD)), 1),'-r','LineWidth',2);
    hold on;
    plot(sf, mean((20*log10(2*noMotionAvgPSD + 0.00001)), 1),'-b','LineWidth',2);
    ylabel('power (dB)');
    xlabel('spatial freq. (cpd)');
    set(gca,'fontsize',16,'yscale','linear','xscale','log');
    xlim([.3 60])
    grid on;
    
    
    subplot(1,7,7)
    plot(tf, mean((20*log10(2*avgPSD)), 2),'-r','LineWidth',2);
    hold on;
    plot(tf, mean((20*log10(2*noMotionAvgPSD + 0.00001)), 2),'-b','LineWidth',2);
    ylabel('power (dB)');
    xlabel('temporal freq. (Hz)');
    set(gca,'fontsize',16,'yscale','linear','xscale','log');
    xlim([2 100])
    grid on;

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