function tempFunc

close all force;
clc;

load('E:\\Eye movement data for Whitening study\\processedData.mat');
minLength = 1.7; % seconds
samplewindow = 128;
pixelSizeDeg = 1/20; 

group = [data.group]; %#ok<NODEF>
dataLength = zeros(length(data),1);
for i=1:length(data)
    dataLength(i) = max(data(i).stitchedTime);
end
usefulData = find(dataLength > minLength);

driftF = [];
driftPS = [];
positionF = [];
positionPs = [];
trialCount = 1;
for i=1:length(usefulData)
    currentTrial = data(usefulData(i));
    Fs = round(1/median(diff(currentTrial.stitchedTime)));
    if Fs~=480
        currentTrial = ResampleTraces(currentTrial);
        data(usefulData(i)) = currentTrial;
    end
    [~,driftF,driftPS(:,trialCount)] =...
        ComputeFFT(480,currentTrial.stitchedPosition, samplewindow);
    [~,positionF,positionPS(:,trialCount)] =...
        ComputeFFT(480,currentTrial.filteredPosition, samplewindow);
    trialCount = trialCount + 1;
end

group = group(usefulData);
groupLabels = {'Old - Drifts','Young - Drifts','Patient - Drifts','Old - Drifts+Microsaccades','Young - Drifts+Microsaccades','Patients - Drifts+Microsaccades'};
colors = [0.93 0.69 0.13
          0 0.45 0.74
          1 0 0.4];
figure;
for i=1:3
    indices = group == i;
    dPS = mean(driftPS(:,indices),2);
    pPS = mean(positionPS(:,indices),2);
    d(i) = semilogx((driftF),10*log10(dPS),'-','Color',colors(i,:),'LineWidth',2); hold on;
    p(i) = semilogx((positionF),10*log10(pPS),'-.','Color',d(i).Color,'LineWidth',2); hold on;
end

xlabel('Temporal frequency (Hz)')
ylabel('Spectral density (dB)')
set(gca,'FontSize',14);
grid on;
legend([d p],groupLabels,'FontSize',10);

% figure;
% edges = 5:1:30;
% for i=1:3
%     indices = [data.group] == i;
%     histogram(dataLength(indices),edges);
%     hold on;
% end
% legend('Old','Young','Patients')
% title(sprintf('Minimum length of trial: %d seconds',minLength));


directory = 'E:\\Eye movement data for Whitening study\\Natural Images database\\To be analyzed\\';
listing = dir(directory);
imageSize = 512;

% go over images
powerSpectra = [];
powerSpectra2 = [];
imageCounter = 1;
try
wb = waitbar(0,'Please wait, images are being analyzed...');
% fh = figure(123);cla;
% set(fh,'units','normalized','outerposition',[.1 .6 .9 .4]);

% this will be needed many many times. so do it here. do it once.
indices = GetIndicesOfMasks(imageSize,imageSize,(samplewindow/2)+1);

st = tic;
for i=19:length(listing)
    if ~(listing(i).isdir)
        [~,~,ext] = fileparts(listing(i).name);
        if ~isempty(strfind(ext,'jpg'))
            
            % read the image
            fullfilename = [directory listing(i).name];
            im = imread(fullfilename);
%             im = imresize(im,4);
            
            % compute the spectral content of the image first
            [m,n]=size(im);
            apertureSizeDeg = imageSize*40/sqrt(sum(size(im).^2));
            startIndex = [(m/2-imageSize/2-1) (n/2-imageSize/2-1)];
            croppedIm = im(startIndex(1):startIndex(1)+imageSize-1, startIndex(2):startIndex(2)+imageSize-1);
            temp = fftshift(fft2(double(croppedIm)));
            [powerSpectra(imageCounter,:), sf] =...
                GetAverageAcrossOrientation((abs(temp).^2)/(imageSize^2), apertureSizeDeg); %#ok<*AGROW>
            
%             figure;
%             semilogx(sf,10*log10(ps)); hold on;
            
            [powerSpectra2(imageCounter,:), sf] =...
                GetAverageAcrossOrientation2((abs(temp).^2)/(imageSize^2), apertureSizeDeg, indices);
            
            
%             semilogx(sf,10*log10(ps))
%             xlabel('Spatial frequency (cpd)')
%             ylabel('Spectral density (db)')
%             set(gca,'FontSize',14);
%             grid on;
%             legend('lines','circles')
            
%             figure(fh);
%             subplot(1,3,1)
%             semilogx(sf,10*log10(powerSpectra(imageCounter,:)),'-b','Linewidth',1);
%             hold on;
%             xlabel('Spatial frequency (cpd)')
%             ylabel('Spectral density ')
%             set(gca,'FontSize',14);
%             grid on;
            
%             % now make movies with eye movements and the current image and
%             % compute 3D-FFT, take average across spatial orientations and
%             % save the resultant 2D spectral density.
%             
% %             poolobj = gcp('nocreate'); % If no pool, do not create new one.
% %             if isempty(poolobj)
% %                 poolobj = parpool(8);
% %             end
%             
%             for j=1:length(usefulData)
%                 currentTrial = data(usefulData(j));
%                 newPixelSizeDeg = apertureSizeDeg/imageSize;
%                 driftShifts = round(currentTrial.stitchedPosition/newPixelSizeDeg);
%                 positionShifts = round(currentTrial.filteredPosition/newPixelSizeDeg);
%                 
%                 % first drifts only
%                 [drift2DPS, driftSF, driftTF, driftFlag] = ...
%                     MakeMovieAndComputePSD(driftShifts, im, imageSize,samplewindow, apertureSizeDeg, indices);
% % %                 figure(fh);
% % %                 subplot(1,3,2);
% % %                 dsuccess = Plot2DPSD(drift2DPS, driftSF, driftTF);
% % %                 title('Drifts only')
%                 
%                 % now all position array: drifts + microsaccades
%                 [position2dPS, posSF, posTF, posFlag] = ...
%                     MakeMovieAndComputePSD(positionShifts, im, imageSize,samplewindow, apertureSizeDeg, indices);
% % %                 figure(fh);
% % %                 subplot(1,3,3);
% % %                 psuccess = Plot2DPSD(position2dPS, posSF, posTF);
% % %                 title('Drifts + Miscrosaccades')
%                 
%                 % create the analysis folder to save files, if it does not
%                 % already exist.
%                 savepath = [fileparts(directory) 'Spectral Analysis New'];
%                 if ~exist(savepath,'dir')
%                     mkdir(savepath);
%                 end
%                     
%                 % save only when both PSD computations completed without a
%                 % problem.
%                 savefilename = [savepath filesep sprintf('im_%d_group_%d_j_%d_noresize_T128.mat',...
%                     imageCounter,currentTrial.group,j)];
%                 currentGroup = currentTrial.group;
%                 trialNumber = j; %#ok<*NASGU>
%                 
%                 SaveFunction(savefilename,drift2DPS,position2dPS,...
%                     driftSF,driftTF,driftFlag,posSF,posTF,posFlag,currentGroup,...
%                     imageCounter,trialNumber);
%                 
%                 fprintf('\n----------------2D PSD saved.-----------------------\n\n');
%             end
            
            imageCounter = imageCounter + 1;
            
            
        end
    end
    waitbar(i/length(listing),wb,sprintf('Progress: %.1f. ETA: %d seconds',...
        i*100/length(listing),round(toc(st)*(length(listing)-i)/i)));
end
catch imerr
    imerr.stack.name
    imerr.stack.line
    delete(wb);
end

try
delete(wb);
catch
end

% try
%     lastf = figure;
%     semilogx(sf,mean(10*log10(powerSpectra),1),'-b','Linewidth',2);
%     hold on;
%     semilogx(sf,mean(10*log10(powerSpectra2),1),'-r','Linewidth',2);
%     semilogx(sf,prctile(10*log10(powerSpectra),2.5,1),'-b');
%     semilogx(sf,prctile(10*log10(powerSpectra),97.5,1),'-b');
%     semilogx(sf,prctile(10*log10(powerSpectra2),2.5,1),'-r');
%     semilogx(sf,prctile(10*log10(powerSpectra2),97.5,1),'-r');
%     xlabel('Spatial frequency (cpd)')
%     ylabel('Spectral density (dB)')
%     set(gca,'FontSize',14);
%     grid on;
% catch
%     close(lastf);
% end
% 
psImagesFilename = [directory 'PSimages_512.mat']; 
save(psImagesFilename,'powerSpectra','powerSpectra2','sf','apertureSizeDeg');

try
delete(poolobj);
catch
end



function SaveFunction(savefilename,drift2DPS,position2dPS,...
    driftSF,driftTF,driftFlag,posSF,posTF,posFlag,currentGroup,...
    imageCounter,trialNumber) %#ok<*INUSD>

save(savefilename,'drift2DPS','position2dPS',...
    'driftSF','driftTF','driftFlag','posSF','posTF','posFlag','currentGroup',...
    'imageCounter','trialNumber');


function success = Plot2DPSD(ps2D, sf, tf)

try
    [ssf, ttf] = meshgrid(sf,tf);
    logPSD = 10*log10(ps2D);
    surf(ssf,ttf,logPSD,'EdgeColor','none');
    view(2)
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    xlim([sf(2) sf(end)]);
    ylim([tf(2) tf(end)])
    colormap(hot);
    xlabel('Spatial frequency (cpd)')
    ylabel('Temporal frequency (Hz)')
    set(gca,'FontSize',14)
    c = colorbar;
    c.Label.String = 'Spectral density (dB)';
    c.FontSize = 10;
    caxis([min(logPSD(:)) max(logPSD(:))]);
    success = 1;
catch ploterr;
    ploterr.stack.name
    ploterr.stack.line
    success = 0;
end


function [ps2D, sf, tf, flag] = MakeMovieAndComputePSD(shifts, im, imageSize, windowSize, apertureSizeDeg, indices)
% this function computes the 2D power spectral density of retinal input
% movies.
% it uses a sliding window of windowSize and overlap specified by increment.
% e.g., if increment = windowSize/4, it means 75% overlap. 
% e.g., if increment = windowSize/2, the overlap is 50%
% NOTE that windowSize is the temporal window size in terms of number of
% samples and imageSize is the spatial window size.

increment = round(windowSize/1);
L = windowSize; % modified to reflect number of samples rather than duration
if rem(L,2)~=0
    L = L-1;
end
howManyTimes = length(1:increment:(length(shifts)-L));
if howManyTimes == 0
    howManyTimes = 1;
end
if howManyTimes > 10;
    howManyTimes = 10;
end

% this light median filtering is needed to reduce the quantization errors
% due to conversion to pixels from visual degrees.
shifts(:,1) = medfilt1(shifts(:,1),3);
shifts(:,2) = medfilt1(shifts(:,2),3);

% create the arrays in advance
movArray = zeros(imageSize,imageSize,windowSize,'uint8','gpuArray');
ps3D = zeros(imageSize,imageSize,windowSize);
im = gpuArray((im));
startIndex = size(im)/2;

% get size info
[m,n,k] = size(movArray);
cx = round(n/2)+1;
cy = round(m/2)+1;

% create necessary arrays
ps2D = zeros((k/2)+1,cx-1);

sf = (0:cx-2)/apertureSizeDeg;
tf = 480*(0:(k/2))/k;

% wb = waitbar(0,'Creating the ps2D array...');
% set(wb,'units','normalized','outerposition',[.4 .2 .2 .1])

accumPSD = zeros((windowSize/2)+1,imageSize/2);
s0 = tic;
try
    for iter = 1 : howManyTimes
%         waitbar(iter/howManyTimes,wb,'Movie is being created...');
        st = tic;
        startTemporalIndex = (iter-1)*increment+1;

        for i=startTemporalIndex:startTemporalIndex+windowSize-1
            movArray(:,:,i-startTemporalIndex+1) = (im(startIndex+shifts(i,2):startIndex+shifts(i,2)+imageSize-1,...
                startIndex+shifts(i,1):startIndex+shifts(i,1)+imageSize-1));
        end
        fprintf('Movie created in %.1f seconds.\n',toc(st));

%         waitbar(iter/howManyTimes,wb,'Now it''s time to compute 2D PSD.');
        
        % do the fft
        st = tic;
        ps3D = gather(abs(ifftshift(fftshift(fftn(double(movArray))),3)).^2)/(m*n*k);
        ps3D = ps3D(:,:,1:(k/2+1));
        ps3D(:,:,2:end-1) = 2*ps3D(:,:,2:end-1);
        fprintf('FFT of the movie took %.1f seconds.\n',toc(st));

        % create empty 2d PS array
        ps2D = zeros(size(ps3D,3),cx-1);
        st = tic;
        
        % average across orientations
        for j=1:size(ps3D,3)
            currentFrame = (ps3D(:,:,j));
            ps2D(j,:) = GetAverageAcrossOrientation2(currentFrame, apertureSizeDeg, indices);   
%             ps2D(j,:) = GetAverageAcrossOrientation(currentFrame, apertureSizeDeg);
        end
        fprintf('Averaging across orientations took %.1f seconds.\n',toc(st));

        accumPSD = accumPSD + ps2D;
    end
    flag = [0 howManyTimes];
    
    clear movArray ps3D
    
    fprintf('\n\n\n\n-----------------TRIAL DONE in %d seconds------------------\n\n\n',round(toc(s0)));
    
catch accumerr
    accumerr.message
    accumerr.stack.name
    accumerr.stack.line
    fprintf('*********Error during 2d PSD*********\niter:%d\n********\n',iter);
    howManyTimes = iter-1;
    flag = [1 howManyTimes];
end

ps2D = accumPSD/howManyTimes;

% delete(wb);





function [ps,f] = GetAverageAcrossOrientation(fftim, apertureSizeDeg)

[m,n] = size(fftim);
cx = round(n/2)+1;
cy = round(m/2)+1;

% figure;

psAll = nan(cx-1,360);
for i=1:360
    [x,y] = pol2cart(i*pi/180,1:cx-2);
    x = [0 x];
    y = [0 y];
    x = round(x) + cx ;
    y = round(y) + cy ;
    
%     plot(x,y); hold on;
    try
        linearInd = sub2ind([m, n], y, x);
    catch inderr
        inderr.stack.name
        inderr.stack.line
        error('GetAverageAcrossOrientation failed.');
    end
    psAll(:,i) = fftim(linearInd)';
end

ps = mean(psAll,2);
f = (0:cx-2)/apertureSizeDeg;

% figure;
% loglog(f,ps)
% xlabel('Spatial frequency (cpd)')
% ylabel('Spectral density ')
% set(gca,'FontSize',14);
% grid on;


function [ps,f] = GetAverageAcrossOrientation2(fftim, apertureSizeDeg, indices)

cx = round(size(fftim,2)/2)+1;
ps = zeros(cx-1,1);
for i=1:(cx-1)
    ps(i) = mean(mean(fftim(indices{i})))';
end

f = (0:cx-2)/apertureSizeDeg;


function indices = GetIndicesOfMasks(m,n,k)

cx = round(n/2)+1;
cy = round(m/2)+1;

inc = m*n;
[y,x] = meshgrid(1:m,1:n);
distance = ((x-cx).^2+(y-cy).^2);
indices = cell(cx-1,1);

% a = zeros(m,n);
for i=1:(cx-1) 
    indices{i} = find((distance >= (i-1)^2) & (distance < i^2));
%     a(indices{i})=rem(i,2)*255;
end

% figure;
% imshow(uint8(a));



function ps2D = GetAverageAcrossOrientation3(fftMov, apertureSizeDeg)

[m,n,k] = size(fftMov);
cx = round(n/2)+1;
cy = round(m/2)+1;

if existsOnGPU(fftMov)
    [y,x] = meshgrid(gpuArray(1:size(fftMov,1)),gpuArray(1:size(fftMov,2)));
    ps2D = zeros(k,cx-1,'gpuArray');
else 
    [y,x] = meshgrid(1:size(fftMov,1),1:size(fftMov,2));
    ps2D = zeros(k,cx-1);
end


for i=1:(cx-1)
    if existsOnGPU(fftMov)
        mask = false([m n],'gpuArray');
    else
        mask = false([m n]);
    end
    distance = ((x-cx).^2+(y-cy).^2);
    indices = (distance >= (i-1)^2) & (distance < i^2);
    mask(indices) = true;

    ps2D(:,i) = squeeze(mean(mean(fftMov.*repmat(mask,1,1,k),1),2))';
%     ps2D(:,i) = squeeze(mean(mean(pagefun(@times,mask,fftMov),1),2))';
    
end



if existsOnGPU(fftMov)
    ps2D = gather(ps2D);
end





function currentTrial = ResampleTraces(currentTrial)

newRawTime = currentTrial.rawTime(1):1/480:currentTrial.rawTime(end);
newStitchedTime = currentTrial.rawTime(1):1/480:currentTrial.rawTime(end);

newRawPosition(:,1) = interp1(currentTrial.rawTime,...
    currentTrial.rawPosition(:,1),newRawTime,'pchip');
newRawPosition(:,2) = interp1(currentTrial.rawTime,...
    currentTrial.rawPosition(:,2),newRawTime,'pchip');
newStitchedPosition(:,1) = interp1(currentTrial.stitchedTime,...
    currentTrial.stitchedPosition(:,1),newStitchedTime,'pchip');
newStitchedPosition(:,2) = interp1(currentTrial.stitchedTime,...
    currentTrial.stitchedPosition(:,2),newStitchedTime,'pchip');

currentTrial.stitchedPosition = newStitchedPosition;
currentTrial.stitchedTime = newStitchedTime;
currentTrial.rawPosition = newRawPosition;
currentTrial.rawTime = newRawTime;


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
