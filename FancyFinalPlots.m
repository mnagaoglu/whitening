function FancyFinalPlots
% The final function to plot everything in publication/presentation quality.

close all force;
clc;

addpath([pwd filesep 'utils/'])
addpath([pwd filesep 'data/'])

backgroundColor = [1 1 1];
axesColor = [0 0 0];
myColors = {'r','k','b'};
% myColors = {'g','w','y'};
alp = 0.6; 
% fileFormat = '-depsc';
fileFormat = '-dtiff';
global fsize;
fsize = 15;

isSave = 0;
if isSave
    if exist([pwd filesep 'results/stats_new.txt'],'file')
        delete([pwd filesep 'results/stats_new.txt']);
    end
    diary([pwd filesep 'results/stats_new.txt']);
end


load('PSimages_512.mat');
try
    powerSpectra = powerSpectra2; % circle averaging method
catch
end

% the following variables will be loaded: 
%       'driftAVG','driftSE','posAVG','posSE',...
%       'numOfObservers','repetitions','drift_psd','pos_psd','sf','tf',...
%       'edriftF','epositionF','edPS','epPS','edPSse','epPSse',...
%       'drift_1d_AVG','drift_1d_SE','pos_1d_AVG','pos_1d_SE','drift_1d_AVG_tf',...
%       'drift_1d_SE_tf','pos_1d_AVG_tf','pos_1d_SE_tf'
load('FinalResultsClean.mat')
% load('intermediateResults.mat')

% get patient PRL data
initials = {'beh';'daj';'dam';'den';'egw';...
            'jrl';'lep';'tgf';'anh';'drm';...
            'eks';'msa';'pfm';'rjm';'jmj'};
PRL = [12.07;8.41;12.73;2.70;1.39;...
       4.10 ;3.57;1.32 ;2.01;5.36;...
       1.66 ;3.39;4.68 ;1.77;10.33]; % deg 
visualAcuity = [1.16; 1.02; 1.1;  0.5;  0.54;
                0.82; 0.62; 0.48; 0.32; 0.98;
                0.74; 0.68; 0.58; 0.5;  0.94];

% reorder PRLs so that values correspond to the list under subjects array.
for i=1:length(subjects{3}) %#ok<*USENS>
    for j=1:length(initials)
        if strcmpi(subjects{3}(i), initials(j))
            newPRL(i) = PRL(j);
            newVisualAcuity(i) = visualAcuity(j);
        end
    end
end
PRL = newPRL;
visualAcuity = newVisualAcuity;

% % FS order is based on subjects{3} array
% fixStab(:,3) = [1.008; 5.377; 2.927; 4.49;0.75;0.324;0.088;2.4;6.49;1.57;8.71;3.84; 2.41;0.57;0.37];


%% 2D plots
%% drifts only
sfl = length(sf);
tfl = length(tf);
figure('units','normalized','outerposition',[.2 .2 .6 .35],...
    'Color',backgroundColor,'InvertHardcopy','off');
subplot(1,3,2);
plot3dInPerspective(sf(2:sfl),tf(2:tfl),...
    10*log10(driftAVG((2:tfl),(2:sfl),2)),0,1,0,'Age-matched',myColors{2});

ax = gca;
ax.Color = backgroundColor;
ax.XColor = axesColor;
ax.YColor = axesColor;

subplot(1,3,1);
plot3dInPerspective(sf(2:sfl),tf(2:tfl),...
    10*log10(driftAVG((2:tfl),(2:sfl),1)),0,0,1,'Young',myColors{2});
ax = gca;
ax.Color = backgroundColor;
ax.XColor = axesColor;
ax.YColor = axesColor;
subplot(1,3,3);
plot3dInPerspective(sf(2:sfl),tf(2:tfl),...
    10*log10(driftAVG((2:tfl),(2:sfl),3)),1,0,0,'AMD',myColors{2});
ax = gca;
ax.Color = backgroundColor;
ax.XColor = axesColor;
ax.YColor = axesColor;
if isSave 
print([pwd filesep 'results' filesep 'drift_2D_new'],'-r300',fileFormat);
end

%% position
figure('units','normalized','outerposition',[.2 .2 .6 .35],...
    'Color',backgroundColor,'InvertHardcopy','off');
subplot(1,3,2);
plot3dInPerspective(sf(2:sfl),tf(2:tfl),...
    10*log10(posAVG((2:tfl),(2:sfl),2)),0,1,0,'Age-matched',myColors{2});
ax = gca;
ax.Color = backgroundColor;
ax.XColor = axesColor;
ax.YColor = axesColor;
subplot(1,3,1);
plot3dInPerspective(sf(2:sfl),tf(2:tfl),...
    10*log10(posAVG((2:tfl),(2:sfl),1)),0,0,1,'Young',myColors{2});
ax = gca;
ax.Color = backgroundColor;
ax.XColor = axesColor;
ax.YColor = axesColor;
subplot(1,3,3);
plot3dInPerspective(sf(2:sfl),tf(2:tfl),...
    10*log10(posAVG((2:tfl),(2:sfl),3)),1,0,0,'AMD',myColors{2});
ax = gca;
ax.Color = backgroundColor;
ax.XColor = axesColor;
ax.YColor = axesColor;
if isSave
print([pwd filesep 'results' filesep 'drift_microsac_2D_new'],'-r300',fileFormat);
end

%% 1D plots
%% temporal
figure('units','normalized','outerposition',[.2 .2 .38 .38],...
    'Color',backgroundColor,'InvertHardcopy','off');
subplot(1,2,1)
th = 1.5;
PlotPatchedSE(tf(2:tfl),(drift_1d_AVG_tf(:,2)),(drift_1d_SE_tf(:,2)),myColors{2},alp); %#ok<*NODEF>
semilogx(tf(2:tfl),10*log10(drift_1d_AVG_tf(:,2)),'-','Color',myColors{2},'LineWidth',th)
hold on;
PlotPatchedSE(tf(2:tfl),(drift_1d_AVG_tf(:,1)),(drift_1d_SE_tf(:,1)),myColors{1},alp);
semilogx(tf(2:tfl),10*log10(drift_1d_AVG_tf(:,1)),'-','Color',myColors{1},'LineWidth',th)
PlotPatchedSE(tf(2:tfl),(drift_1d_AVG_tf(:,3)),(drift_1d_SE_tf(:,3)),myColors{3},alp);
semilogx(tf(2:tfl),10*log10(drift_1d_AVG_tf(:,3)),'-','Color',myColors{3},'LineWidth',th)
xlabel('Temporal frequency (hz)')
ylabel('Spectral density (dB)')
title('drifts only')
xlim([4 240])
ylim([0 60])
grid on;
set(gca,'FontSize',fsize)
set(gca,'XTick',[10 100]);
set(gca,'XTickLabel',{'10','100'})
ax = gca;
ax.Color = backgroundColor;
ax.XColor = axesColor;
ax.YColor = axesColor;


subplot(1,2,2)
PlotPatchedSE(tf(2:tfl),(pos_1d_AVG_tf(:,2)),(pos_1d_SE_tf(:,2)),myColors{2},alp);
a(1) = semilogx(tf(2:tfl),10*log10(pos_1d_AVG_tf(:,2)),'-','Color',myColors{2},'LineWidth',th);
hold on;
PlotPatchedSE(tf(2:tfl),(pos_1d_AVG_tf(:,1)),(pos_1d_SE_tf(:,1)),myColors{1},alp);
a(2) = semilogx(tf(2:tfl),10*log10(pos_1d_AVG_tf(:,1)),'-','Color',myColors{1},'LineWidth',th);
PlotPatchedSE(tf(2:tfl),(pos_1d_AVG_tf(:,3)),(pos_1d_SE_tf(:,3)),myColors{3},alp);
a(3) = semilogx(tf(2:tfl),10*log10(pos_1d_AVG_tf(:,3)),'-','Color',myColors{3},'LineWidth',th);
xlabel('Temporal frequency (hz)')
% ylabel('Spectral density (dB)')
lh= legend(a,{'Young','Age-Matched','AMD'},'Location','best');
lh.FontSize = 10;
lh.Color = backgroundColor;
lh.TextColor = myColors{2};
lh.EdgeColor = myColors{2};
xlim([4 240])
ylim([0 60])
grid on;
title('drifts + microsaccades')
set(gca,'FontSize',fsize)
set(gca,'XTick',[10 100]);
set(gca,'XTickLabel',{'10','100'})
ax = gca;
ax.Color = backgroundColor;
ax.XColor = axesColor;
ax.YColor = axesColor;

if isSave
print([pwd filesep 'results' filesep 'temporal_1D_new'],'-r300',fileFormat);
end




%% spatial 
figure('units','normalized','outerposition',[.2 .2 .38 .38],...
    'Color',backgroundColor,'InvertHardcopy','off');
subplot(1,2,1)
PlotPatchedSE(sf(2:sfl),(drift_1d_AVG(:,2)),(drift_1d_SE(:,2)),myColors{2},alp);
semilogx(sf(2:sfl),10*log10(drift_1d_AVG(:,2)),'-','Color',myColors{2},'LineWidth',th)
hold on;
PlotPatchedSE(sf(2:sfl),(drift_1d_AVG(:,1)),(drift_1d_SE(:,1)),myColors{1},alp);
semilogx(sf(2:sfl),10*log10(drift_1d_AVG(:,1)),'-','Color',myColors{1},'LineWidth',th)
PlotPatchedSE(sf(2:sfl),(drift_1d_AVG(:,3)),(drift_1d_SE(:,3)),myColors{3},alp);
semilogx(sf(2:sfl),10*log10(drift_1d_AVG(:,3)),'-','Color',myColors{3},'LineWidth',th)
xlabel('Spatial Frequency (cpd)')
ylabel('Spectral density (dB)')
grid on;
semilogx(sf(2:sfl),3+10*log10(mean(powerSpectra(:,2:end),1)),'--','Color',myColors{2},'LineWidth',th)
xlim([.2 65])
ylim([0 70])
title('drifts only')
set(gca,'FontSize',fsize)
set(gca,'XTick',[1 10])
set(gca,'XTickLabel',{'1','10'})
ax = gca;
ax.Color = backgroundColor;
ax.XColor = axesColor;
ax.YColor = axesColor;

subplot(1,2,2)
PlotPatchedSE(sf(2:sfl),(pos_1d_AVG(:,2)),(pos_1d_SE(:,2)),myColors{2},alp);
a(1) = semilogx(sf(2:sfl),10*log10(pos_1d_AVG(:,2)),'-','Color',myColors{2},'LineWidth',th);
hold on;
PlotPatchedSE(sf(2:sfl),(pos_1d_AVG(:,1)),(pos_1d_SE(:,1)),myColors{1},alp);
a(2) = semilogx(sf(2:sfl),10*log10(pos_1d_AVG(:,1)),'-','Color',myColors{1},'LineWidth',th);
PlotPatchedSE(sf(2:sfl),(pos_1d_AVG(:,3)),(pos_1d_SE(:,3)),myColors{3},alp);
a(3) = semilogx(sf(2:sfl),10*log10(pos_1d_AVG(:,3)),'-','Color',myColors{3},'LineWidth',th);
xlabel('Spatial Frequency (cpd)')
% ylabel('Spectral density (dB)')
a(4) = semilogx(sf(2:sfl),3+10*log10(mean(powerSpectra(:,2:end),1)),'--','Color',myColors{2},'LineWidth',th);
lh = legend(a,{'Young','Age-Matched','AMD','Images'},'Location','best');
lh.FontSize = 10;
lh.Color = backgroundColor;
lh.TextColor = myColors{2};
lh.EdgeColor = myColors{2};
grid on;
xlim([.2 65])
ylim([0 70])
title('drifts + microsaccades')
set(gca,'FontSize',fsize)
set(gca,'XTick',[1 10])
set(gca,'XTickLabel',{'1','10'})
ax = gca;
ax.Color = backgroundColor;
ax.XColor = axesColor;
ax.YColor = axesColor;


if isSave
print([pwd filesep 'results' filesep 'spatial_1D_new'],'-r300',fileFormat);
end


%% natural scene power coefficient
means = sum(powerSpectra(:,2:end),2);
[~,I] = sort(means,'descend');
powerSpectra = powerSpectra(I,:);
cols = jet(size(powerSpectra,1));

if ~exist('imagePS.mat','file')
    for i=1:size(powerSpectra,1)
        maxi = find(sf>10,1);
        currentPS = log10(powerSpectra(i,2:maxi));
        currentSF = (sf(2:maxi));
        imCoeff(i) = getPowerCoeff(log10(currentSF), currentPS);
    end
    save('imagePS.mat','imCoeff','currentSF','maxi');
else
    load('imagePS.mat');
end


figure('units','normalized','outerposition',[.2 .2 .20 .38],...
    'Color',backgroundColor,'InvertHardcopy','off');
histogram(imCoeff,15,'FaceColor',[.6 .6 .6],'EdgeColor',myColors{2});
xlabel('Power Coefficient')
ylabel('Frequency')
box off
set(gca,'FontSize',fsize)
xlim([0.5 3.5])
ylim([0 15])
ax = gca;
ax.Color = backgroundColor;
ax.XColor = axesColor;
ax.YColor = axesColor;

if isSave
    print([pwd filesep 'results' filesep 'imagePower_new'],'-r300',fileFormat);
end


%% compute and plot whitening factors
for i=1:length(subjects)
    numOfObservers(i) = length(subjects{i});
end

for j=1:3
    for i=1:numOfObservers(j)
        driftCoeff(i,j) = getPowerCoeff(log10(currentSF), log10(drift_1d_psd_sf(2:maxi,i,j))');
        posCoeff(i,j) = getPowerCoeff(log10(currentSF), log10(pos_1d_psd_sf(2:maxi,i,j))');
    end
end




figure('units','normalized','outerposition',[.2 .2 .38 .38],...
    'Color',backgroundColor,'InvertHardcopy','off');
msize = 40;
m=0.25;

whitened = mean(imCoeff);
youngWF = (whitened-driftCoeff(1:numOfObservers(2),2))/whitened;
agedWF = (whitened-driftCoeff(1:numOfObservers(1),1))/whitened;
patientWF = (whitened-driftCoeff(1:numOfObservers(3),3))/whitened;
youngWF_pos = (whitened-posCoeff(1:numOfObservers(2),2))/whitened;
agedWF_pos = (whitened-posCoeff(1:numOfObservers(1),1))/whitened;
patientWF_pos = (whitened-posCoeff(1:numOfObservers(3),3))/whitened;

fprintf('\n-----------------------------------------------------------\n');
fprintf('\n\n\n Two-sample t-tests for drifts only\n\n')
[h_PatientVsAged,p_PatientVsAged,ci_PatientVsAged,stats_PatientVsAged] = ...
    ttest2(patientWF,agedWF)
fprintf('\n-----------------------------------------------------------\n');
[h_youngVsAged,p_youngVsAge,ci_youngVsAge,stats_youngVsAge] = ...
    ttest2(youngWF,agedWF)
fprintf('\n-----------------------------------------------------------\n');
[h_youngVsPatient,p_youngVsPatient,ci_youngVsPatient,stats_youngVsPatient] = ...
    ttest2(youngWF,patientWF)
fprintf('\n-----------------------------------------------------------\n');
fprintf('\n Two-sample t-tests for drifts + microsaccades\n\n')
[h_PatientVsAged,p_PatientVsAged,ci_PatientVsAged,stats_PatientVsAged] = ...
    ttest2(patientWF_pos,agedWF_pos)
fprintf('\n-----------------------------------------------------------\n');
[h_youngVsAged,p_youngVsAge,ci_youngVsAge,stats_youngVsAge] = ...
    ttest2(youngWF_pos,agedWF_pos)
fprintf('\n-----------------------------------------------------------\n');
[h_youngVsPatient,p_youngVsPatient,ci_youngVsPatient,stats_youngVsPatient] = ...
    ttest2(youngWF_pos,patientWF_pos)
fprintf('\n-----------------------------------------------------------\n');

seY = std(youngWF)/sqrt(sum(~isnan(youngWF)));
seA = std(agedWF)/sqrt(sum(~isnan(agedWF)));
seP = std(patientWF)/sqrt(sum(~isnan(patientWF)));
seY_pos = std(youngWF_pos)/sqrt(sum(~isnan(youngWF_pos)));
seA_pos = std(agedWF_pos)/sqrt(sum(~isnan(agedWF_pos)));
seP_pos = std(patientWF_pos)/sqrt(sum(~isnan(patientWF_pos)));
patchAlpha = 0.1;
pha = patch([mean(youngWF)-seY mean(youngWF)+seY mean(youngWF)+seY mean(youngWF)-seY],...
    [0 0 12 12],myColors{2},'FaceColor',myColors{2},'EdgeColor','none'); hold on;
alpha(pha,patchAlpha);
pha = patch([mean(agedWF)-seA mean(agedWF)+seA mean(agedWF)+seA mean(agedWF)-seA],...
    [0 0 12 12],myColors{1},'FaceColor',myColors{1},'EdgeColor','none'); hold on;
alpha(pha,patchAlpha);
pha = patch([mean(patientWF)-seP mean(patientWF)+seP mean(patientWF)+seP mean(patientWF)-seP],...
    [0 0 12 12],myColors{3},'FaceColor',myColors{3},'EdgeColor','none'); hold on;
alpha(pha,patchAlpha);
pha = patch([mean(youngWF_pos)-seY_pos mean(youngWF_pos)+seY_pos mean(youngWF_pos)+seY_pos mean(youngWF_pos)-seY_pos],...
    [0 0 12 12],myColors{2},'FaceColor',myColors{2},'EdgeColor','none'); hold on;
alpha(pha,patchAlpha);
pha = patch([mean(agedWF_pos)-seA_pos mean(agedWF_pos)+seA_pos mean(agedWF_pos)+seA_pos mean(agedWF_pos)-seA_pos],...
    [0 0 12 12],myColors{1},'FaceColor',myColors{1},'EdgeColor','none'); hold on;
alpha(pha,patchAlpha);
pha = patch([mean(patientWF_pos)-seP_pos mean(patientWF_pos)+seP_pos mean(patientWF_pos)+seP_pos mean(patientWF_pos)-seP_pos],...
    [0 0 12 12],myColors{3},'FaceColor',myColors{3},'EdgeColor','none'); hold on;
alpha(pha,patchAlpha);


plot(mean(youngWF)*[1 1],[0 12],'--','Color',myColors{2},'LineWidth',1);hold on;
plot(mean(agedWF)*[1 1],[0 12],'--','Color',myColors{1},'LineWidth',1);
plot(mean(patientWF)*[1 1],[0 12],'--','Color',myColors{3},'LineWidth',1);
plot(mean(youngWF_pos)*[1 1],[0 12],'-','Color',myColors{2},'LineWidth',1);
plot(mean(agedWF_pos)*[1 1],[0 12],'-','Color',myColors{1},'LineWidth',1);
plot(mean(patientWF_pos)*[1 1],[0 12],'-','Color',myColors{3},'LineWidth',1);

y = randn(numOfObservers(2),1)*m;
hh(1) = scatter((whitened-driftCoeff(1:numOfObservers(2),2))/whitened,y+3,msize,myColors{2},'MarkerEdgeColor',myColors{2},'Marker','o'); hold on;
% alpha(hh(1),alp/1000);
y = randn(numOfObservers(1),1)*m;
sh = scatter((whitened-driftCoeff(1:numOfObservers(1),1))/whitened,y+7,msize,myColors{1},'MarkerEdgeColor',myColors{1},'Marker','s');
% alpha(sh,alp/1000);
y = randn(numOfObservers(3),1)*m;
sh =scatter((whitened-driftCoeff(1:numOfObservers(3),3))/whitened,y+11,msize,myColors{3},'MarkerEdgeColor',myColors{3},'Marker','d');
% alpha(sh,alp/1000);

y = randn(numOfObservers(2),1)*m;
hh(2) = scatter((whitened-posCoeff(1:numOfObservers(2),2))/whitened,y+3,msize,myColors{2},'filled','MarkerEdgeColor',myColors{2},'Marker','o'); hold on;
alpha(hh(1),alp);
y = randn(numOfObservers(1),1)*m;
sh = scatter((whitened-posCoeff(1:numOfObservers(1),1))/whitened,y+7,msize,myColors{1},'filled','MarkerEdgeColor',myColors{1},'Marker','s');
alpha(sh,alp);
y = randn(numOfObservers(3),1)*m;
sh = scatter((whitened-posCoeff(1:numOfObservers(3),3))/whitened,y+11,msize,myColors{3},'filled','MarkerEdgeColor',myColors{3},'Marker','d');
alpha(sh,alp);

lh = legend(hh,{'drifts only','drifts + microsaccades'},'Location','best');
lh.FontSize = 10;
lh.Color = backgroundColor;
lh.TextColor = myColors{2};
lh.EdgeColor = myColors{2};

xlabel('"Whitening" factor')
set(gca,'FontSize',fsize)
set(gca,'YTick',[3 7 11])
set(gca,'YTickLabel',{'Young','Age-Matched','AMD'})
xlim([0 1.2])
ylim([0 16])
box off
ax = gca;
ax.Color = backgroundColor;
ax.XColor = axesColor;
ax.YColor = axesColor;

if isSave
print([pwd filesep 'results' filesep 'whiteningFactors_new'],'-r300',fileFormat);
end



%% eye movement spectra
%    'edriftF','epositionF','edPS','epPS','edPSse','epPSse'

figure('units','normalized','outerposition',[.2 .2 .38 .38],...
    'Color',backgroundColor,'InvertHardcopy','off');
subplot(1,2,1)
th = 1.5;
xl = length(edriftF);
PlotPatchedSE(edriftF(2:xl),(edPS(2:xl,2)),(edPSse(2:xl,2)),myColors{2},alp);
semilogx(edriftF(2:xl),10*log10(edPS(2:xl,2)),'-','Color',myColors{2},'LineWidth',th)
hold on;
PlotPatchedSE(edriftF(2:xl),(edPS(2:xl,1)),(edPSse(2:xl,1)),myColors{1},alp);
semilogx(edriftF(2:xl),10*log10(edPS(2:xl,1)),'-','Color',myColors{1},'LineWidth',th)
PlotPatchedSE(edriftF(2:xl),(edPS(2:xl,3)),(edPSse(2:xl,3)),myColors{3},alp);
semilogx(edriftF(2:xl),10*log10(edPS(2:xl,3)),'-','Color',myColors{3},'LineWidth',th)
xlabel('Temporal frequency (hz)')
ylabel('Spectral density (dB)')
title('drifts only')
xlim([4 240])
ylim([-60 20])
grid on;
set(gca,'FontSize',fsize)
set(gca,'XTick',[3 10 30 100]);
set(gca,'XTickLabel',{'3','10','30','100'})
ax = gca;
ax.Color = backgroundColor;
ax.XColor = axesColor;
ax.YColor = axesColor;

subplot(1,2,2)
PlotPatchedSE(edriftF(2:xl),(epPS(2:xl,2)),(epPSse(2:xl,2)),myColors{2},alp);
a(1) = semilogx(edriftF(2:xl),10*log10(epPS(2:xl,2)),'-','Color',myColors{2},'LineWidth',th);
hold on;
PlotPatchedSE(edriftF(2:xl),(epPS(2:xl,1)),(epPSse(2:xl,1)),myColors{1},alp);
a(2) = semilogx(edriftF(2:xl),10*log10(epPS(2:xl,1)),'-','Color',myColors{1},'LineWidth',th);
PlotPatchedSE(edriftF(2:xl),(epPS(2:xl,3)),(epPSse(2:xl,3)),myColors{3},alp);
a(3) = semilogx(edriftF(2:xl),10*log10(epPS(2:xl,3)),'-','Color',myColors{3},'LineWidth',th);
xlabel('Temporal frequency (hz)')
% ylabel('Spectral density (dB)')
lh= legend(a(1:3),{'Young','Age-Matched','AMD'});
lh.FontSize = 10;
lh.Color = backgroundColor;
lh.TextColor = myColors{2};
lh.EdgeColor = myColors{2};
xlim([4 240])
ylim([-60 20])
grid on;
title('drifts + miscrosaccades')
set(gca,'FontSize',fsize)
set(gca,'XTick',[3 10 30 100]);
set(gca,'XTickLabel',{'3','10','30','100'})
ax = gca;
ax.Color = backgroundColor;
ax.XColor = axesColor;
ax.YColor = axesColor;


if isSave
print([pwd filesep 'results' filesep 'eyemovementSpectra_new'],'-r300',fileFormat);
end



%% PRL vs WF
figure('Color',backgroundColor,'InvertHardcopy','off','units','normalized',...
    'outerposition',[.1 .1 0.1990    0.3861]);
% subplot(1,4,1)
patientW = (whitened-driftCoeff(1:numOfObservers(3),3))/whitened;
% scatter(PRL, patientW); hold on;
tbl = table(PRL',patientW,'VariableNames',{'PRL','patientW'});
lm = fitlm(tbl,'patientW~PRL')
ph = lm.plot;
ph(1).MarkerSize = 6;
ph(1).Marker = 'd';
ph(1).MarkerEdgeColor = myColors{3};
ph(1).MarkerFaceColor = myColors{3};
alpha(alp);
ph(2).Color = myColors{3};
ph(2).LineWidth = 2;
ph(3).Color = myColors{3};
ph(3).LineWidth = 2;
ph(4).Color = myColors{3};
ph(4).LineWidth = 2;
set(gca,'FontSize',fsize)
box off;
xlabel(sprintf('PRL eccentricity \n (deg)'))
ylabel('Whitening factor')
title('')
lh = findobj(gcf,'Tag','legend');
lh.Location = 'best';
lh.FontSize = 9;
lh.String{3} = 'CI';
lh.Color = backgroundColor;
lh.TextColor = myColors{2};
lh.EdgeColor = myColors{2};
ax = gca;
ax.Color = backgroundColor;
ax.XColor = axesColor;
ax.YColor = axesColor;

if isSave
print([pwd filesep 'results' filesep 'correlations_PRL'],'-r300',fileFormat);
end


%% Fixation stability vs Whitening factor
% subplot(1,4,3)
figure('Color',backgroundColor,'InvertHardcopy','off','units','normalized',...
    'outerposition',[.1 .1 0.1990    0.3861]);
youngW = (whitened-driftCoeff(1:numOfObservers(2),2))/whitened;
oldW = (whitened-driftCoeff(1:numOfObservers(1),1))/whitened;
patientW = (whitened-driftCoeff(1:numOfObservers(3),3))/whitened;

fsAll = [fixStab(1:numOfObservers(2),2); fixStab(1:numOfObservers(1),1); fixStab(1:numOfObservers(3),3)];
wfAll = [youngW; oldW; patientW];

tbl2 = table(log10(fsAll),wfAll,'VariableNames',{'FixationStability','WhiteningFactor'});
lm2 = fitlm(tbl2,'WhiteningFactor~FixationStability')
tblPatient = table(log10(fixStab(1:numOfObservers(3),3)),patientW,'VariableNames',{'FixationStability','WhiteningFactor'});
lmPatient = fitlm(tblPatient,'WhiteningFactor~FixationStability')

sh = scatter(log10(fixStab(1:numOfObservers(2),2)), youngW,msize,myColors{2},...
    'filled','MarkerEdgeColor',myColors{2},'Marker','o'); hold on;
alpha(sh,alp);
sh = scatter(log10(fixStab(1:numOfObservers(1),1)), oldW,msize,myColors{1},...
    'filled','MarkerEdgeColor',myColors{1},'Marker','s');
alpha(sh,alp);
sh = scatter(log10(fixStab(1:numOfObservers(3),3)), patientW,msize,myColors{3},...
    'filled','MarkerEdgeColor',myColors{3},'Marker','d');
alpha(sh,alp);
plot(log10(linspace(min(fsAll),max(fsAll),100)),...
    lm2.feval(log10(linspace(min(fsAll),max(fsAll),100))),...
    '-','Color',[.4 .4 .4],'LineWidth',2);
plot(log10(linspace(min(fixStab(1:numOfObservers(3),3)),max(fixStab(1:numOfObservers(3),3)),100)),...
    lmPatient.feval(log10(linspace(min(fixStab(1:numOfObservers(3),3)),max(fixStab(1:numOfObservers(3),3)),100))),...
    '-','Color',myColors{3},'LineWidth',2);
set(gca,'FontSize',fsize)
box off;
xlabel(sprintf('Fixation instability \n (log[deg^2])'))
ylabel('Whitening factor')
lh = legend('Young','Age-Matched','AMD','Fit to all','Fit to AMD');
lh.Location = 'best';
lh.FontSize = 9;
lh.Color = backgroundColor;
lh.TextColor = myColors{2};
lh.EdgeColor = myColors{2};
% xlim([-0.5 3.5])
ax = gca;
ax.Color = backgroundColor;
ax.XColor = axesColor;
ax.YColor = axesColor;

if isSave
print([pwd filesep 'results' filesep 'correlations_fixstab'],'-r300',fileFormat);
end



%% Visual acuity vs whitening factor
figure('Color',backgroundColor,'InvertHardcopy','off','units','normalized',...
    'outerposition',[.1 .1 0.1990    0.3861]);
% subplot(1,4,2)
patientW = (whitened-driftCoeff(1:numOfObservers(3),3))/whitened;
tbl = table(visualAcuity',patientW,'VariableNames',{'visualAcuity','patientW'});
lm3 = fitlm(tbl,'patientW~visualAcuity')
ph = lm3.plot;
ph(1).MarkerSize = 6;
ph(1).Marker = 'd';
ph(1).MarkerEdgeColor = myColors{3};
ph(1).MarkerFaceColor = myColors{3};
alpha(alp)
ph(2).Color = myColors{3};
ph(2).LineWidth = 2;
ph(3).Color = myColors{3};
ph(3).LineWidth = 2;
ph(4).Color = myColors{3};
ph(4).LineWidth = 2;
set(gca,'FontSize',fsize)
box off;
xlabel(sprintf('Visual acuity \n (logMAR)'))
ylabel('Whitening factor')
title('')
lh = findobj(gcf,'Tag','legend');
lh = lh(1);
lh.Location = 'best';
lh.FontSize = 9;
lh.String{3} = 'CI';
lh.Color = backgroundColor;
lh.TextColor = myColors{2};
lh.EdgeColor = myColors{2};
ax = gca;
ax.Color = backgroundColor;
ax.XColor = axesColor;
ax.YColor = axesColor;

if isSave
print([pwd filesep 'results' filesep 'correlations_va'],'-r300',fileFormat);
end


%% Diffusion vs WF
load(['../ganglion/diffusionData_n128.mat']);
% reorder diffusion constants so that values correspond to the list under subjects array.
for i=1:length(subjects{3}) %#ok<*USENS>
    for j=1:length(diffusionData)
        if strcmpi(subjects{3}(i), diffusionData(j).initials) && ...
                (diffusionData(j).group == 3)
            diffusionDrift(i) = diffusionData(j).D_fit_drift;
            diffusionPos(i) = diffusionData(j).D_fit;
        end
    end
end

allSubjects = [subjects{2}'; subjects{1}'; subjects{3}'];
allW = [youngW; oldW; patientW];
groupNumbers = [2*ones(size(youngW)); 1*ones(size(oldW)); 3*ones(size(patientW))];
for i=1:length(allSubjects) %#ok<*USENS>
    for j=1:length(diffusionData)
        if strcmpi(allSubjects{i}, diffusionData(j).initials) && ...
                (groupNumbers(i) == diffusionData(j).group)
            allDrift(i) = diffusionData(j).D_fit_drift;
            allPos(i) = diffusionData(j).D_fit;
        end
    end
end



figure('Color',backgroundColor,'InvertHardcopy','off','units','normalized',...
    'outerposition',[.1 .1 0.1990    0.3861]);
% subplot(1,4,4);
cla;
tblAllDiffusion = table(log10(allDrift)',allW,'VariableNames',{'Diffusion','WhiteningFactor'});
lmAllDiffusion = fitlm(tblAllDiffusion,'WhiteningFactor~Diffusion')
tblDiffPatient = table(log10(diffusionDrift'),patientW,'VariableNames',{'Diffusion','WhiteningFactor'});
lmDiffPatient = fitlm(tblDiffPatient,'WhiteningFactor~Diffusion')


sh = scatter(log10(allDrift(1:numOfObservers(2))), youngW,msize,myColors{2},...
    'filled','MarkerEdgeColor',myColors{2},'Marker','o'); hold on;
alpha(sh,alp);
sh = scatter(log10(allDrift(numOfObservers(2)+1:sum(numOfObservers(1:2)))), oldW,msize,myColors{1},...
    'filled','MarkerEdgeColor',myColors{1},'Marker','s');
alpha(sh,alp);
sh = scatter(log10(allDrift(sum(numOfObservers(1:2))+1:end)), patientW,msize,myColors{3},...
    'filled','MarkerEdgeColor',myColors{3},'Marker','d');
alpha(sh,alp);

plot(log10(linspace(min(allDrift),max(allDrift),100)),...
    lmAllDiffusion.feval(log10(linspace(min(allDrift),max(allDrift),100))),...
    '-','Color',[.4 .4 .4],'LineWidth',2);
plot(log10(linspace(min(diffusionDrift),max(diffusionDrift),100)),...
    lmDiffPatient.feval(log10(linspace(min(diffusionDrift),max(diffusionDrift),100))),...
    '-','Color',myColors{3},'LineWidth',2);
set(gca,'FontSize',fsize)
box off;
xlabel(sprintf('log[Diffusion constant] \n (log[arcmin^2/sec])'))
ylabel('Whitening factor')
lh = legend('Young','Age-Matched','AMD','Fit to all','Fit to AMD');
lh.Location = 'best';
lh.FontSize = 9;
lh.Color = backgroundColor;
lh.TextColor = myColors{2};
lh.EdgeColor = myColors{2};
% xlim([-0.5 3.5])
ax = gca;
ax.Color = backgroundColor;
ax.XColor = axesColor;
ax.YColor = axesColor;

if isSave
print([pwd filesep 'results' filesep 'correlations_diffusion'],'-r300',fileFormat);
end













if isSave
    diary off;
end

end



function PlotPatchedSE(x,y,se,c,a)

if isrow(y)
    y = y';
end
if isrow(se)
    se = se';
end
if isrow(x)
    x = x';
end
py = [y+se; flipud(y-se)];
px = [x; flipud(x)];

pp = patch(px,10*log10(py),c);
pp.FaceAlpha = a;
pp.EdgeColor = 'none';
set(gca,'XScale','log')
hold on;
end



function plot3dInPerspective(sf,tf,S, ...
    isColorbar,isX,isY,titlestr,cols)

if nargin<4
    isColorbar = 0;
end

newSf = logspace(min(log10(sf)),max(log10(sf)),512);
newTf = logspace(min(log10(tf)),max(log10(tf)),512);
[X,Y] = meshgrid(sf,tf);
[Xq,Yq] = meshgrid(newSf,newTf);
newS = interp2(X,Y,S,Xq,Yq,'cubic');


surf(newSf,newTf,newS,'EdgeColor','none','LineWidth',0.25,'EdgeAlpha',.5);
xlim([.25 65])
ylim([4 240])
minZ = -5;
maxZ = 72;
zlim([minZ maxZ])
view(2);set(gca,'XScale','log');
set(gca,'YScale','log');
colormap(jet(256))
set(gca,'CLim',[minZ maxZ])
set(gca,'XTick',[1 3 10 30])
set(gca,'XTickLabel',{'1','3','10','30'})
set(gca,'YTick',[3 10 30 100]);
set(gca,'YTickLabel',{'3','10','30','100'})
ax = gca;
ax.TickDir = 'out';
ax.TickLength = [0.02 0.02];
title(titlestr);

if isX
    xlabel('Spatial frequency (cpd)');
else
    xlabel(' ')
end
if isY
	ylabel('Temporal frequency (Hz)');
else
    ylabel(' ')
end
global fsize;
set(gca,'FontSize',fsize);
pause(.2);
s1Pos = get(gca,'position');
if isColorbar
    cb = colorbar('eastoutside');
    cb.Label.String = 'Spectral Density (dB)';
    cb.Color = cols;
end
pause(.2);
s2Pos = get(gca,'position');
s2Pos(3:4) = s1Pos(3:4);
set(gca,'position',s2Pos);

if isColorbar
   cb.Position = cb.Position./[1 1 3 1]; 
end


end



function [coeff, fitresult, gof] = getPowerCoeff(currentSF, currentPS)
%CREATEFIT(CURRENTSF,CURRENTPS)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : currentSF
%      Y Output: currentPS
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 01-Feb-2017 16:34:17


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( currentSF, currentPS );

% Set up fittype and options.
ft = fittype( '(-a*x+b)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Algorithm = 'Levenberg-Marquardt';
opts.Display = 'Off';
opts.Robust = 'Bisquare';
opts.StartPoint = [0.722905352239109 0.811753323578174];


% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

coeff = fitresult.a;
% 
% % Plot fit with data.
% figure(99);hold on;
% h = plot( fitresult, xData, yData );
% legend( h, 'currentPS vs. currentSF', 'untitled fit 1', 'Location', 'NorthEast' );
% % Label axes
% xlabel currentSF
% ylabel currentPS
% grid on
% set(gca,'XScale','log')
% set(gca,'YScale','log')


end