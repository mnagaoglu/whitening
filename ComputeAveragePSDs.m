function ComputeAveragePSDs
close all
clc

% mypath = 'E:\\Eye movement data for Whitening study\\Natural Images database\\Spectral Analysis\\';
% mypath = 'E:\Eye movement data for Whitening study\Natural Images database\To be analyzed\Spectral Analysis New\';
mypath = 'E:\Eye movement data for Whitening study\Natural Images database\To be analyzed\Spectral Analysis With Clean EM\';

% load('E:\Eye movement data for Whitening study\processedData.mat')
% missing = load('E:\Eye movement data for Whitening study\missingData.mat');
% 
% data = [data; missing.data];
% clear missing;
% load('FS.mat');

load('C:\Users\spencer\Google Drive\SELAB code\Mehmet\Whitening Analysis\ganglion\readyToAnalyzeData_new.mat')
data = eyeMovements;
clear eyeMovements;

% FS = FixStab(data);
load('FS_new.mat','FS');

%% eye movement spectra
minLength = 1.7; % seconds
samplewindow = 2048;

group = [data.group]; 
dataLength = zeros(length(data),1);
for i=1:length(data)
    dataLength(i) = max(data(i).stitchedTime);
end

% % remove the following data entries
% toremove = [27 41 51 62 64 67 75:80 89:91 94:96 99:100 114:115 128:130 132];
toInclude = true(length(data),1);
% toInclude(toremove) = false;

usefulData = find((dataLength > minLength) & toInclude);

edriftF = [];
edriftPS = [];
epositionF = [];
epositionPS = [];
trialCount = 1;
for i=1:length(usefulData)
    currentTrial = data(usefulData(i));
    Fs = round(1/median(diff(currentTrial.stitchedTime)));
    if Fs~=480
        currentTrial = ResampleTraces(currentTrial);
        data(usefulData(i)) = currentTrial;
    end
    [~,edriftF,edriftPS(:,trialCount)] =...
        ComputeFFT(480,currentTrial.stitchedPosition, samplewindow);
    [~,epositionF,epositionPS(:,trialCount)] =...
        ComputeFFT(480,currentTrial.filteredPosition, samplewindow);
    trialCount = trialCount + 1;
end

group = group(usefulData);

for i=1:3
    indices = group == i;
    edPS(:,i) = mean(edriftPS(:,indices),2);
    epPS(:,i) = mean(epositionPS(:,indices),2);
    edPSse(:,:,i) = std(edriftPS(:,indices),[],2)/sqrt(sum(indices));
    epPSse(:,:,i) = std(epositionPS(:,indices),[],2)/sqrt(sum(indices));
end



%% retinal power spectra
initials = [data.initials];
initials = upper(reshape(initials,3,length(data))');
for i=1:length(initials)
    init{i} = initials(i,:);
end

groups = [data.group];
numOfObservers = [];
for i=1:3
    indices = (groups == i) & toInclude';
    numOfObservers(i) = length(unique(init(indices)));
    subjects{i} = unique(init(indices));
end

T = 65;
N = 256;
drift_psd = zeros(T,N,max(numOfObservers),3);
pos_psd = zeros(T,N,max(numOfObservers),3);
repetitions = zeros(max(numOfObservers),3);
fixStab = repetitions;

drift_1d_psd_sf = zeros(N-1,max(numOfObservers),3);
pos_1d_psd_sf = zeros(N-1,max(numOfObservers),3);

drift_1d_psd_tf = zeros(T-1,max(numOfObservers),3);
pos_1d_psd_tf = zeros(T-1,max(numOfObservers),3);


clear data;


listing = dir(mypath);

% figure;
cols = [1 0 0; 0 0 1; 0 0 0];
wb = waitbar(0,'Please wait, files are being loaded..');
for i=1:length(listing)
    if ~(listing(i).isdir)
        [~,~,ext] = fileparts(listing(i).name);
        if ~isempty(strfind(ext,'mat'))
            fullfilename = [mypath listing(i).name];
            load(fullfilename);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             trialNumber = trialNumber-125;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            try
            [observerIndex,groupIndex,goAhead] = ...
                FindIndices(init,groups,trialNumber,toInclude);
            
            if goAhead
                drift_psd(:,:,observerIndex,groupIndex) = drift_psd(:,:,observerIndex,groupIndex) + drift2DPS;
                pos_psd(:,:,observerIndex,groupIndex) = pos_psd(:,:,observerIndex,groupIndex) + position2dPS;
                repetitions(observerIndex,groupIndex) = repetitions(observerIndex,groupIndex) + 1;
                fixStab(observerIndex,groupIndex) =  FS(trialNumber);
                sf = driftSF;
                tf = driftTF;


                % average across nonzero temporal or spatial frequencies
                drift_1d_psd_sf(:,observerIndex,groupIndex) = drift_1d_psd_sf(:,observerIndex,groupIndex) + mean(drift2DPS(2:end,2:end),1)';
                pos_1d_psd_sf(:,observerIndex,groupIndex) = pos_1d_psd_sf(:,observerIndex,groupIndex) + mean(position2dPS(2:end,2:end),1)';
                drift_1d_psd_tf(:,observerIndex,groupIndex) = drift_1d_psd_tf(:,observerIndex,groupIndex) + mean(drift2DPS(2:end,2:end),2);
                pos_1d_psd_tf(:,observerIndex,groupIndex) = pos_1d_psd_tf(:,observerIndex,groupIndex) + mean(position2dPS(2:end,2:end),2);
    %             semilogx(sf(2:end),10*log10(spatialPS(2:end)),'-','Color',cols(groupIndex,:));hold on;
            end
            catch errr
                errr.message
                errr.stack.name
                errr.stack.line
            end
        end
    end
    waitbar(i/length(listing),wb,sprintf('Progress %d',round(100*i/length(listing))))
end

% xlabel('Spatial frequency (cpd)')
% ylabel('Power spectra (dB)')
% xlim([.15 65]);grid on;
% set(gca,'FontSize',14) 
        
delete(wb);

%% average across repetititons and images
for i =1:size(repetitions,1)
    for j=1:size(repetitions,2)
        drift_psd(:,:,i,j) = drift_psd(:,:,i,j)/repetitions(i,j);
        pos_psd(:,:,i,j) = pos_psd(:,:,i,j)/repetitions(i,j); 
        drift_1d_psd_sf(:,i,j) = drift_1d_psd_sf(:,i,j)/repetitions(i,j);
        pos_1d_psd_sf(:,i,j) =  pos_1d_psd_sf(:,i,j)/repetitions(i,j);
        drift_1d_psd_tf(:,i,j) =  drift_1d_psd_tf(:,i,j)/repetitions(i,j);
        pos_1d_psd_tf(:,i,j) =  pos_1d_psd_tf(:,i,j)/repetitions(i,j);
    end
end

%% average across observers
driftAVG = zeros(T,N,3);
driftSE = driftAVG;
posAVG = zeros(T,N,3);
posSE = posAVG;

drift_1d_AVG = zeros(N-1,3);
drift_1d_SE = drift_1d_AVG;
pos_1d_AVG = drift_1d_AVG;
pos_1d_SE = drift_1d_AVG;
drift_1d_AVG_tf = zeros(T-1,3);
drift_1d_SE_tf = drift_1d_AVG_tf;
pos_1d_AVG_tf = drift_1d_AVG_tf;
pos_1d_SE_tf = drift_1d_AVG_tf;


for i=1:3
    driftAVG(:,:,i) = squeeze(mean(drift_psd(:,:,1:numOfObservers(i),i),3));
    driftSE(:,:,i) = squeeze(std(drift_psd(:,:,1:numOfObservers(i),i),0,3))/sqrt(numOfObservers(i));
    posAVG(:,:,i) = squeeze(mean(pos_psd(:,:,1:numOfObservers(i),i),3));
    posSE(:,:,i) = squeeze(std(pos_psd(:,:,1:numOfObservers(i),i),0,3))/sqrt(numOfObservers(i));
    
    drift_1d_AVG(:,i) = squeeze(mean(drift_1d_psd_sf(:,1:numOfObservers(i),i),2));
    drift_1d_SE(:,i) = squeeze(std(drift_1d_psd_sf(:,1:numOfObservers(i),i),0,2))/sqrt(numOfObservers(i));
    pos_1d_AVG(:,i) = squeeze(mean(pos_1d_psd_sf(:,1:numOfObservers(i),i),2));
    pos_1d_SE(:,i) = squeeze(std(pos_1d_psd_sf(:,1:numOfObservers(i),i),0,2))/sqrt(numOfObservers(i));
    drift_1d_AVG_tf(:,i) = squeeze(mean(drift_1d_psd_tf(:,1:numOfObservers(i),i),2));
    drift_1d_SE_tf(:,i) = squeeze(std(drift_1d_psd_tf(:,1:numOfObservers(i),i),0,2))/sqrt(numOfObservers(i));
    pos_1d_AVG_tf(:,i) = squeeze(mean(pos_1d_psd_tf(:,1:numOfObservers(i),i),2));
    pos_1d_SE_tf(:,i) = squeeze(std(pos_1d_psd_tf(:,1:numOfObservers(i),i),0,2))/sqrt(numOfObservers(i));
end


save('FinalResultsClean.mat','driftAVG','driftSE','posAVG','posSE',...
    'numOfObservers','repetitions','drift_psd','pos_psd','sf','tf',...
    'edriftF','epositionF','edPS','epPS','edPSse','epPSse',...
    'drift_1d_AVG','drift_1d_SE','pos_1d_AVG','pos_1d_SE','drift_1d_AVG_tf',...
    'drift_1d_SE_tf','pos_1d_AVG_tf','pos_1d_SE_tf',...
    'drift_1d_psd_sf','pos_1d_psd_sf','drift_1d_psd_tf','pos_1d_psd_tf','subjects','fixStab');


        

function [Y,f,PS] = ComputeFFT(Fs,position, minLength)
% this function computes the power spectral density of eye position traces
% it uses a sliding window of minLength and overlap specified by increment.
% e.g., if increment = minLength/4, it means 75% overlap. 
% e.g., if increment = minLength/2, the overlap is 50%


increment = round(minLength/2);
L = minLength; % modified to reflect number of samples rather than duration
if rem(L,2)~=0
    L = L-1;
end
howManyTimes = length(1:increment:(length(position)-L));
if howManyTimes == 0
    howManyTimes = 1;
end

for i=1:howManyTimes
    iterPosition = sqrt(sum(position((i-1)*increment+1:(i-1)*increment+L,:).^2,2));
    Y(:,i) = fft(iterPosition);
    P2 = abs(Y).^2/L;
    P1 = P2(1:(L/2+1));
    P1(2:end-1) = 2*P1(2:end-1);
    PS(:,i) = P1;
end

Y = mean(Y,2);
PS = mean(PS,2);
f = Fs*(0:(L/2))/L;
     

function currentTrial = ResampleTraces(currentTrial)

newRawTime = currentTrial.rawTime(1):1/480:currentTrial.rawTime(end);
newStitchedTime = currentTrial.rawTime(1):1/480:currentTrial.rawTime(end);

newRawPosition(:,1) = interp1(currentTrial.rawTime,...
    currentTrial.rawPosition(:,1),newRawTime,'pchip',NaN);
newRawPosition(:,2) = interp1(currentTrial.rawTime,...
    currentTrial.rawPosition(:,2),newRawTime,'pchip',NaN);
newStitchedPosition(:,1) = interp1(currentTrial.stitchedTime,...
    currentTrial.stitchedPosition(:,1),newStitchedTime,'pchip',NaN);
newStitchedPosition(:,2) = interp1(currentTrial.stitchedTime,...
    currentTrial.stitchedPosition(:,2),newStitchedTime,'pchip',NaN);

currentTrial.stitchedPosition = newStitchedPosition;
currentTrial.stitchedTime = newStitchedTime;
currentTrial.rawPosition = newRawPosition;
currentTrial.rawTime = newRawTime;


function [o,g, goAhead] = FindIndices(initials, groups, trialNumber, toInclude)
goAhead = toInclude(trialNumber);

if goAhead
    g = groups(trialNumber);
    observers = unique(initials((g == groups) & toInclude'));
    thisObserver = initials(trialNumber);
    for i=1:length(observers)
        if strcmp(observers(i),thisObserver)
            o=i;
            return;
        end
    end
else
    o = [];
    g = [];
end



function FS = FixStab(data)

for i=1:length(data)
    pos = data(i).filteredPosition;
    indices = data(i).rawTime < 30;
    FS(i) = GetISOA(pos(indices,:));
end


function ISOA = GetISOA(pos)

[density, X, Y, bandwidth, PRL, fh, stats] = GetKSDensity(pos(:,1), pos(:,2), 0);

ISOA = stats.areaInContour;



