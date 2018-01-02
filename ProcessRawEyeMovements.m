function ProcessRawEyeMovements
% Process eye movement data for further processing. needed only once.

clearvars;
close all force;

toSaveData = 'E:\\Eye movement data for Whitening study\\processedData.mat';

folders{1} = 'E:\\Eye movement data for Whitening study\\Analysed Data\\Normals (Old)\\';
folders{2} = 'E:\\Eye movement data for Whitening study\\Analysed Data\\Normals (Young)\\';
folders{3} = 'E:\\Eye movement data for Whitening study\\Analysed Data\\Patients (UCB)\\';
folders{4} = 'E:\\Eye movement data for Whitening study\\Analysed Data\\Patients (UH)\\';

% some parameters
keyPhrases = {'480_hz_';'480_hz_';'540_hz_';'480_hz_'};
samplingRate = [480; 480; 540; 480];
pixelSizeDeg = 1/20; 



% walk through the folders listed above 
data = [];
for i=1:length(folders)
    currentListing = dir(folders{i});
    usefulFiles = zeros(size(currentListing));
    
    % refine the list by filtering out the irrelevant files/folders
    for j=1:length(currentListing)
       if ~currentListing(j).isdir
           if ~isempty(strfind(currentListing(j).name,keyPhrases(i)))
               usefulFiles(j) = 1;
            
               % while at it, load the relevant file
               fileContents = load([folders{i} currentListing(j).name]);
               
               % find saccades and plot results, make a selection as to
               % whether or not use this file for further analysis
               [isGood, saccades, rawPosition, rawTime, filteredPosition,...
                   stitchedPosition, stitchedTime]...
                   = PreProcess(fileContents, pixelSizeDeg,samplingRate(i));
               
               if isGood == 1
                   thisTrial.saccades = saccades;
                   thisTrial.rawPosition = rawPosition;
                   thisTrial.rawTime = rawTime;
                   thisTrial.filteredPosition = filteredPosition;
                   thisTrial.stitchedPosition = stitchedPosition;
                   thisTrial.stitchedTime = stitchedTime;
                   filename = currentListing(j).name; 
                   thisTrial.fullFileName = [folders{i} filename];
                   thisTrial.initials = filename(1:3);
                   if i>2
                       group = 3;
                   else
                       group = i;
                   end
                   thisTrial.group = group; 
                   data = [data; thisTrial];
               end
               
               if isGood == -1
                   return;
               end
               
           end
       end
    end
    
end
clear currentListing;

save(toSaveData,'data');



function [isGood, saccades, rawPosition, rawTime, filteredPosition,...
    stitchedPosition, stitchedTime] ...
    = PreProcess(fileContents, pixelSizeDeg, sampleRate)

try
    time = fileContents.timeaxis_secs;
    position = fileContents.frameshifts_strips_spline * pixelSizeDeg;

    % shift everything to zero time and position
    time = time - time(1);
    position(:,1) = position(:,1) - position(1,1);
    position(:,2) = position(:,2) - position(1,2);

    % fix the gaps in time
    [time, position] = FixTemporalGaps(time, position, sampleRate);

    % find out the number of strips per frame and do median filtering with a
    % window of half this size.
    numberOfStripsPerFrame = fileContents.numbervideoframes / fileContents.videoframerate;
    mfposition(:,1) = medfilt1(position(:,1),round(numberOfStripsPerFrame/2));
    mfposition(:,2) = medfilt1(position(:,2),round(numberOfStripsPerFrame/2));

    % lowpass filtering eye positions 
    lpFilt = designfilt('lowpassfir','SampleRate',sampleRate,...
             'PassbandFrequency',120, ...
             'StopbandFrequency',150,'PassbandRipple',0.01, ...
             'StopbandAttenuation',300,'DesignMethod','kaiserwin');
    lpmfposition(:,1) = filtfilt(lpFilt,mfposition(:,1));
    lpmfposition(:,2) = filtfilt(lpFilt,mfposition(:,2));
    
    % notch filter
    d = designfilt('bandstopiir','FilterOrder',2, ...
               'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
               'DesignMethod','butter','SampleRate',sampleRate);
    bsmflpposition(:,1) = filtfilt(d,lpmfposition(:,1));
    bsmflpposition(:,2) = filtfilt(d,lpmfposition(:,2));

%     figure(21234);cla;
%     [Y,f] = ComputeFFT(sampleRate,sqrt(position(:,1).^2 + position(:,2).^2));
%     [Y,f] = ComputeFFT(sampleRate,sqrt(mfposition(:,1).^2 + mfposition(:,2).^2));
%     [Y,f] = ComputeFFT(sampleRate,sqrt(lpmfposition(:,1).^2 + lpmfposition(:,2).^2));
%     [Y,f] = ComputeFFT(sampleRate,sqrt(bsmflpposition(:,1).^2 + bsmflpposition(:,2).^2));
    
    % now find saccades in the 2D velocity space
    saccades = FindSaccades(time, bsmflpposition);

    % get drifts
    [stitchedPosition, stitchedTime] = FindDrifts(time,bsmflpposition,saccades, sampleRate);

    filteredPosition = bsmflpposition;
    rawPosition = position;
    rawTime = time;

    isGood = 1;

    figure(12345);
    subplot(4,1,4)
    cla;
    plot(stitchedTime,stitchedPosition(:,1),'-r',stitchedTime,stitchedPosition(:,2),'-b')
    ylabel('Drift only')
    
    % Construct a questdlg with two options
    choice = questdlg('Do you want to use this trial', ...
    	'Is Good Menu', ...
    	'Yes','No','Exit','No');
    % Handle response
    switch choice
        case 'Yes'
            isGood = 1;
        case 'No'
            isGood = 0;
        case 'Exit'
            isGood = -1;   
        otherwise
            isGood = -1;
    end
    
    if isempty(saccades)
        isGood = 1;       
    end
    
catch preperr

    isGood = 0;
    saccades = []; 
    rawPosition = []; 
    rawTime = [];  
    stitchedPosition = [];  
    stitchedTime = []; 
    filteredPosition = [];
    preperr.stack.name
    preperr.stack.line
end



function [Y,f,P1] = ComputeFFT(Fs,position)

L = length(position); % Length of signal

Y = fft(position);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;
figure(21234);
semilogy(f,P1); hold on;
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')



function [newTime, newPosition] = FixTemporalGaps(time, position, sampleRate)
% this function finds the temporal gaps in the data and
% interpolates/stitches gap region appropriately. If the gap is tiny, e.g.,
% a few samples long, then we interpolate.. If longer than that, we stitch
% different parts.

n = 5; % samples
diffTime = diff(time);
indicesTBI = find((diffTime > ((1/sampleRate)+0.0001)) & (diffTime < (n*(1/sampleRate)+0.0001)));

% start interpolation
newTime = time;
newPosition = position;
offsetDueToInsertion = 0;
for i=1:length(indicesTBI)
    tempTime = time(indicesTBI(i):indicesTBI(i)+1);
    tempPosition = position(indicesTBI(i):indicesTBI(i)+1,:);
    newTempTime = (tempTime(1):(1/sampleRate):tempTime(2)+0.0002)';
    newTempPosition(:,1) = interp1(tempTime,tempPosition(:,1),newTempTime,'pchip');
    newTempPosition(:,2) = interp1(tempTime,tempPosition(:,2),newTempTime,'pchip');
    
    % now inset the interpolated regions back to time and position arrays
    newTime = [ newTime(1:(indicesTBI(i)+offsetDueToInsertion-1)); ...
                newTempTime; ...
                newTime((indicesTBI(i)+offsetDueToInsertion+2):end)];
    newPosition = [ newPosition(1:(indicesTBI(i)+offsetDueToInsertion-1),:); ...
                    newTempPosition; ...
                    newPosition((indicesTBI(i)+offsetDueToInsertion+2):end,:)];
    offsetDueToInsertion = offsetDueToInsertion + length(newTempTime) - 2;
    
end

% now to-be-stitched ones
diffTime = diff(newTime);
indicesTBS = find(diffTime >= (n*(1/sampleRate)+0.0001));

for i=1:length(indicesTBS)
    deltaT = -newTime(indicesTBS(i)+1) + newTime(indicesTBS(i)) + (1/sampleRate);
    newTime(indicesTBS(i)+1:end) = newTime(indicesTBS(i)+1:end) + deltaT;
    deltaPosition = -newPosition(indicesTBS(i)+1,:) + newPosition(indicesTBS(i),:);
    newPosition(indicesTBS(i)+1:end,:) = newPosition(indicesTBS(i)+1:end,:) + repmat(deltaPosition,(length(newTime)-indicesTBS(i)),1);
end


% % double check if we got it right
% N = histcounts(diff(newTime));
% if N == (length(newTime)-1)
%     fprintf('Temporal gaps have been fixed successfully!\n');
% else
%     fprintf('There is something fishy with the time array. Check it out more carefully.\n')
% end


function saccades = FindSaccades(time, position, lambda, deltaStitch, nSample)

if nargin<3
    lambda = 6;
    deltaStitch = 10;
    nSample = 8;
end
maxDuration = 0.15; % secs
regionLength = max(time)/2; % seconds

try

    % compute velocity
    velocity_hor = [0; diff(position(:,1))./diff(time)];
    velocity_ver = [0; diff(position(:,2))./diff(time)];
    
    % now detect saccades (from Englbert & Kliegl 2003)
    
    threshold_hor = GetRegionalThreshold(velocity_hor,time,regionLength,lambda);
    threshold_ver = GetRegionalThreshold(velocity_ver,time,regionLength,lambda);
    aboveThreshold = (abs(velocity_hor) > threshold_hor) | (abs(velocity_ver) > threshold_ver);

    % take the difference of indices computed above to find the onset and
    % offset of the movement
    dabove = [0; diff(aboveThreshold)];
    onsets = find(dabove == 1);
    offsets = find(dabove == -1);

    if length(onsets) > length(offsets)
        offsets = [offsets; length(velocity_hor)];
    elseif length(onsets) < length(offsets)
        offsets = offsets(1:end-1);
    end


    % if two consecutive saccades are closer than 10ms, merge them
    tempOnsets = onsets;
    tempOffsets = offsets;

    % loop until there is no pair of consecutive saccades closer than
    % deltaStitch ms
    while true
        for c=1:min(length(onsets),length(offsets))-1
            if (onsets(c+1)-offsets(c))<deltaStitch 
                tempOnsets(c+1) = -1;
                tempOffsets(c) = -1;
            end
        end

        s_on = tempOnsets(tempOnsets ~= -1);
        s_off = tempOffsets(tempOffsets ~= -1);
        if sum((s_on(2:end)-s_off(1:end-1)) < deltaStitch)==0
            break;
        end
    end

    onsets = tempOnsets(tempOnsets ~= -1);
    offsets = tempOffsets(tempOffsets ~= -1);
    
    % remove too brief and too long saccades
    tooBrief = (offsets - onsets) < nSample;
    tooLong = (time(offsets) - time(onsets)) > maxDuration;
    toBeDiscarded = [onsets(tooBrief | tooLong) offsets(tooBrief | tooLong)];
    onsets(tooBrief | tooLong) = [];
    offsets(tooBrief | tooLong) = [];
    
    % duration of saccades
    durations = time(offsets) - time(onsets);

    % finally compute the amplitude of saccades
    amplitudes = abs(position(offsets) - position(onsets));

    saccades.onsets = onsets;
    saccades.offsets = offsets;
    saccades.amplitudes = amplitudes;
    saccades.durations = durations;
    saccades.toBeDiscarded = toBeDiscarded;
    saccades.lambda = lambda;
    saccades.deltaStitch = deltaStitch;
    saccades.nSample = nSample;
    
    
    h=figure(12345);
    set(h,'units','normalized','outerposition',[.05 .1 .5 .8]);
    subplot(4,1,1)
    cla;
    plot(velocity_hor,velocity_ver,'-k');
    hold on;
    for i=1:length(onsets)
        plot(velocity_hor(onsets(i):offsets(i)),velocity_ver(onsets(i):offsets(i)),'-r');
    end
    title('2D velocity space')
    subplot(4,1,2)
    cla;
    plot(time,velocity_hor,'-r',time,velocity_ver,'-b');hold on;
    for i=1:length(onsets)
        plot(time(onsets(i):offsets(i)),velocity_hor(onsets(i):offsets(i)),'or',time(onsets(i):offsets(i)),velocity_ver(onsets(i):offsets(i)),'ob');
    end
    ylabel('Velocity')
    subplot(4,1,3)
    cla;
    plot(time,position(:,1),'-r',time,position(:,2),'-b');hold on;
    for i=1:length(onsets)
        plot(time(onsets(i):offsets(i)),position(onsets(i):offsets(i),1),'or',time(onsets(i):offsets(i)),position(onsets(i):offsets(i),2),'ob');
    end
    ylabel('Position')
    
catch findsacerr

    saccades = [];
    findsacerr.stack.name
    findsacerr.stack.line
end



function threshold = GetRegionalThreshold(velocity,time,regionLength,lambda)


regions = 0:regionLength:max(time);
if regions(end) ~= max(time)
    regions(end) = max(time);
end

threshold = zeros(size(velocity));
for i=1:(length(regions)-1)
    indices = (regions(i)<=time) & (regions(i+1)>=time);
    threshold(indices) = lambda*sqrt(median(velocity(indices).^2)-median(velocity(indices))^2);
end



function [stitchedPosition, stitchedTime] = FindDrifts(time,position,saccades,sampleRate)

% first remove the saccade regions from time and position arrays
onsets = saccades.onsets;
offsets = saccades.offsets;
for i=1:length(onsets)
    time(onsets(i):offsets(i)) = NaN;
    position(onsets(i):offsets(i),:) = NaN;
end

% also discard the toBeDiscarded
for i=1:length(saccades.toBeDiscarded)
    time(saccades.toBeDiscarded(i,1):saccades.toBeDiscarded(i,2)) = NaN;
    position(saccades.toBeDiscarded(i,1):saccades.toBeDiscarded(i,2),:) = NaN;
end
nanIndices = isnan(time);
time(nanIndices) = [];
position(nanIndices,:) = [];

% stitch across saccadic gaps and discarded regions to get drift only traces.. adjust the time
% array accordingly, as if nothing happened.
diffTime = diff(time);
indicesTBS = find(diffTime >= (1/sampleRate)+0.0001);

% now it's time to stitch across saccades
stitchedTime = time;
stitchedPosition = position;
for i=1:length(indicesTBS)
    deltaT = -stitchedTime(indicesTBS(i)+1) + stitchedTime(indicesTBS(i)) + (1/sampleRate);
    stitchedTime(indicesTBS(i)+1:end) = stitchedTime(indicesTBS(i)+1:end) + deltaT;
    deltaPosition = -stitchedPosition(indicesTBS(i)+1,:) + stitchedPosition(indicesTBS(i),:);
    stitchedPosition(indicesTBS(i)+1:end,:) = ...
        stitchedPosition(indicesTBS(i)+1:end,:) + repmat(deltaPosition,(length(stitchedTime)-indicesTBS(i)),1);
end





