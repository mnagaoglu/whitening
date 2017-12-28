function FinalPlots
close all force;
clc;


load('E:\\Eye movement data for Whitening study\\Natural Images database\\To be analyzed\\PSimages_512.mat');

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

% reorder PRLs so that values correspond to the list under subjects array.
for i=1:length(subjects{3}) %#ok<*USENS>
    for j=1:length(initials)
        if strcmpi(subjects{3}(i), initials(j))
            newPRL(i) = PRL(j);
        end
    end
end
PRL = newPRL;




%% 2D plots
%% drifts only
maxval = 60;
minval = 0;
figure;
subplot(1,3,2);
surf(sf(2:end),tf(2:end),10*log10(driftAVG((2:end),(2:end),1)),'EdgeColor','none','LineWidth',0.25,'EdgeAlpha',.5);
xlim([.25 65])
ylim([4 240])
view(2);set(gca,'XScale','log');set(gca,'YScale','log');colormap(jet(500))
set(gca,'CLim',[minval maxval])
set(gca,'XTick',[1 10])
set(gca,'XTickLabel',{'1','10'})
set(gca,'YTick',[10 100]);
set(gca,'YTickLabel',{'10','100'})

xlabel('Spatial frequency (cpd)')
title('Age-matched')
set(gca,'FontSize',14)
subplot(1,3,1);
surf(sf(2:end),tf(2:end),10*log10(driftAVG((2:end),(2:end),2)),'EdgeColor','none','LineWidth',0.25,'EdgeAlpha',.5);
xlim([.25 65])
ylim([4 240])
view(2);set(gca,'XScale','log');set(gca,'YScale','log');colormap(jet(500))
set(gca,'CLim',[minval maxval])
ylabel('Temporal frequency (Hz)')
title('Young')
set(gca,'FontSize',14)
set(gca,'XTick',[1 10])
set(gca,'XTickLabel',{'1','10'})
set(gca,'YTick',[10 100]);
set(gca,'YTickLabel',{'10','100'})

subplot(1,3,3);
surf(sf(2:end),tf(2:end),10*log10(driftAVG((2:end),(2:end),3)),'EdgeColor','none','LineWidth',0.25,'EdgeAlpha',.5);
xlim([.25 65])
ylim([4 240])
view(2);set(gca,'XScale','log');set(gca,'YScale','log');colormap(jet(500))
set(gca,'CLim',[minval maxval])
title('AMD')
set(gca,'FontSize',14)
set(gca,'XTick',[1 10])
set(gca,'XTickLabel',{'1','10'})
set(gca,'YTick',[10 100]);
set(gca,'YTickLabel',{'10','100'})
c=colorbar;
c.Label.String = 'Spectral Density (dB)';
c.Label.FontSize = 12;


%% position

figure;
subplot(1,3,2);
surf(sf(2:end),tf(2:end),10*log10(posAVG((2:end),(2:end),1)),'EdgeColor','none','LineWidth',0.25,'EdgeAlpha',.5);
xlim([.25 65])
ylim([4 240])
view(2);set(gca,'XScale','log');set(gca,'YScale','log');colormap(jet(500))
set(gca,'CLim',[minval maxval])
set(gca,'XTick',[1 10])
set(gca,'XTickLabel',{'1','10'})
set(gca,'YTick',[10 100]);
set(gca,'YTickLabel',{'10','100'})

xlabel('Spatial frequency (cpd)')
title('Age-matched')
set(gca,'FontSize',14)
subplot(1,3,1);
surf(sf(2:end),tf(2:end),10*log10(posAVG((2:end),(2:end),2)),'EdgeColor','none','LineWidth',0.25,'EdgeAlpha',.5);
xlim([.25 65])
ylim([4 240])
view(2);set(gca,'XScale','log');set(gca,'YScale','log');colormap(jet(500))
set(gca,'CLim',[minval maxval])
ylabel('Temporal frequency (Hz)')
title('Young')
set(gca,'FontSize',14)
set(gca,'XTick',[1 10])
set(gca,'XTickLabel',{'1','10'})
set(gca,'YTick',[10 100]);
set(gca,'YTickLabel',{'10','100'})

subplot(1,3,3);
surf(sf(2:end),tf(2:end),10*log10(posAVG((2:end),(2:end),3)),'EdgeColor','none','LineWidth',0.25,'EdgeAlpha',.5);
xlim([.25 65])
ylim([4 240])
view(2);set(gca,'XScale','log');set(gca,'YScale','log');colormap(jet(500))
set(gca,'CLim',[minval maxval])
title('AMD')
set(gca,'FontSize',14)
set(gca,'XTick',[1 10])
set(gca,'XTickLabel',{'1','10'})
set(gca,'YTick',[10 100]);
set(gca,'YTickLabel',{'10','100'})
c=colorbar;
c.Label.String = 'Spectral Density (dB)';
c.Label.FontSize = 12;












%% 1D plots
figure;
subplot(1,2,1)
alp = 0.25;
th = 1.5;
PlotPatchedSE(tf(2:end),(drift_1d_AVG_tf(:,2)),(drift_1d_SE_tf(:,2)),'k',alp);
semilogx(tf(2:end),10*log10(drift_1d_AVG_tf(:,2)),'-k','LineWidth',th)

hold on;
PlotPatchedSE(tf(2:end),(drift_1d_AVG_tf(:,1)),(drift_1d_SE_tf(:,1)),'r',alp);
semilogx(tf(2:end),10*log10(drift_1d_AVG_tf(:,1)),'-r','LineWidth',th)
PlotPatchedSE(tf(2:end),(drift_1d_AVG_tf(:,3)),(drift_1d_SE_tf(:,3)),'b',alp);
semilogx(tf(2:end),10*log10(drift_1d_AVG_tf(:,3)),'-b','LineWidth',th)
xlabel('Temporal frequency (hz)')
ylabel('Spectral density (dB)')
title('drifts only')
xlim([4 240])
grid on;
set(gca,'FontSize',14)

set(gca,'XTick',[10 100]);
set(gca,'XTickLabel',{'10','100'})

subplot(1,2,2)

PlotPatchedSE(tf(2:end),(pos_1d_AVG_tf(:,2)),(pos_1d_SE_tf(:,2)),'k',alp);
a(1) = semilogx(tf(2:end),10*log10(pos_1d_AVG_tf(:,2)),'-k','LineWidth',th);
hold on;

PlotPatchedSE(tf(2:end),(pos_1d_AVG_tf(:,1)),(pos_1d_SE_tf(:,1)),'r',alp);
a(2) = semilogx(tf(2:end),10*log10(pos_1d_AVG_tf(:,1)),'-r','LineWidth',th);

PlotPatchedSE(tf(2:end),(pos_1d_AVG_tf(:,3)),(pos_1d_SE_tf(:,3)),'b',alp);
a(3) = semilogx(tf(2:end),10*log10(pos_1d_AVG_tf(:,3)),'-b','LineWidth',th);
xlabel('Temporal frequency (hz)')
ylabel('Spectral density (dB)')
lh= legend(a,{'Young','Age-Matched','AMD'});
lh.FontSize = 10;
xlim([4 240])
grid on;
title('drifts + microsaccades')
set(gca,'FontSize',14)
set(gca,'XTick',[10 100]);
set(gca,'XTickLabel',{'10','100'})

figure;
subplot(1,2,1)

PlotPatchedSE(sf(2:end),(drift_1d_AVG(:,2)),(drift_1d_SE(:,2)),'k',alp);
semilogx(sf(2:end),10*log10(drift_1d_AVG(:,2)),'-k','LineWidth',th)
hold on;
PlotPatchedSE(sf(2:end),(drift_1d_AVG(:,1)),(drift_1d_SE(:,1)),'r',alp);
semilogx(sf(2:end),10*log10(drift_1d_AVG(:,1)),'-r','LineWidth',th)
PlotPatchedSE(sf(2:end),(drift_1d_AVG(:,3)),(drift_1d_SE(:,3)),'b',alp);
semilogx(sf(2:end),10*log10(drift_1d_AVG(:,3)),'-b','LineWidth',th)
xlabel('Spatial Frequency (cpd)')
ylabel('Spectral density (dB)')
grid on;
semilogx(sf(2:end),3+10*log10(mean(powerSpectra(:,2:end),1)),'--k','LineWidth',th)
xlim([.2 65])
title('drifts only')
set(gca,'FontSize',14)
set(gca,'XTick',[1 10])
set(gca,'XTickLabel',{'1','10'})

subplot(1,2,2)
PlotPatchedSE(sf(2:end),(pos_1d_AVG(:,2)),(pos_1d_SE(:,2)),'k',alp);
a(1) = semilogx(sf(2:end),10*log10(pos_1d_AVG(:,2)),'-k','LineWidth',th);
hold on;
PlotPatchedSE(sf(2:end),(pos_1d_AVG(:,1)),(pos_1d_SE(:,1)),'r',alp);
a(2) = semilogx(sf(2:end),10*log10(pos_1d_AVG(:,1)),'-r','LineWidth',th);
PlotPatchedSE(sf(2:end),(pos_1d_AVG(:,3)),(pos_1d_SE(:,3)),'b',alp);
a(3) = semilogx(sf(2:end),10*log10(pos_1d_AVG(:,3)),'-b','LineWidth',th);
xlabel('Spatial Frequency (cpd)')
ylabel('Spectral density (dB)')
a(4) = semilogx(sf(2:end),3+10*log10(mean(powerSpectra(:,2:end),1)),'--k','LineWidth',th);
lh = legend(a,{'Young','Age-Matched','AMD','Images'});
lh.FontSize = 10;
grid on;
xlim([.2 65])
title('drifts + microsaccades')
set(gca,'FontSize',14)
set(gca,'XTick',[1 10])
set(gca,'XTickLabel',{'1','10'})


%% natural scene power coefficient
% [x,y] = meshgrid(sf,1:length(powerSpectra));

means = sum(powerSpectra(:,2:end),2);
[~,I] = sort(means,'descend');
powerSpectra = powerSpectra(I,:);
cols = jet(size(powerSpectra,1));

% figure(123);
for i=1:size(powerSpectra,1)
    maxi = find(sf>10,1);
    currentPS = log10(powerSpectra(i,2:maxi));
    currentSF = (sf(2:maxi));
    imCoeff(i) = getPowerCoeff(log10(currentSF), currentPS);
%     figure(123);
%     plot(currentSF,10*(currentPS),'Color',cols(i,:));
%     hold on;
end
% xlabel('Spatial frequency')
% ylabel('Power spectra (dB)')
% set(gca,'XScale','log');
% % xlim([.5 100])
% set(gca,'FontSize',14)
% grid on;
numOfObservers = [14,7,15];

for j=1:3
    for i=1:numOfObservers(j)
        driftCoeff(i,j) = getPowerCoeff(log10(currentSF), log10(drift_1d_psd_sf(2:maxi,i,j))');
        posCoeff(i,j) = getPowerCoeff(log10(currentSF), log10(pos_1d_psd_sf(2:maxi,i,j))');
    end
end
figure;
msize = 60;
m=0.25;

whitened = mean(imCoeff);
% y = randn(size(imCoeff))*m/2;
% scatter((whitened-imCoeff)/whitened,y+1,msize/2,'w','filled','MarkerEdgeColor','k'); hold on;
% plot([1 1]*mean((whitened-imCoeff)/whitened),[0 9],'--k','LineWidth',1); hold on; 
y = randn(numOfObservers(2),1)*m;
hh(1) = scatter((whitened-driftCoeff(1:numOfObservers(2),2))/whitened,y+3,msize,[.5 .5 .5],'filled','MarkerEdgeColor','k','Marker','o'); hold on;

y = randn(numOfObservers(1),1)*m;
scatter((whitened-driftCoeff(1:numOfObservers(1),1))/whitened,y+5,msize,[1 .5 .5],'filled','MarkerEdgeColor','r','Marker','s')
y = randn(numOfObservers(3),1)*m;
scatter((whitened-driftCoeff(1:numOfObservers(3),3))/whitened,y+7,msize,[.5 .5 1],'filled','MarkerEdgeColor','b','Marker','d')


y = randn(numOfObservers(2),1)*m;
hh(2) = scatter((whitened-posCoeff(1:numOfObservers(2),2))/whitened,y+3,msize,'w','filled','MarkerEdgeColor','k','Marker','o'); hold on;
y = randn(numOfObservers(1),1)*m;
scatter((whitened-posCoeff(1:numOfObservers(1),1))/whitened,y+5,msize,'w','filled','MarkerEdgeColor','r','Marker','s')
y = randn(numOfObservers(3),1)*m;
scatter((whitened-posCoeff(1:numOfObservers(3),3))/whitened,y+7,msize,'w','filled','MarkerEdgeColor','b','Marker','d')

plot(mean((whitened-driftCoeff(1:numOfObservers(2),2))/whitened)*[1 1],[0 9],'k','LineWidth',2);
plot(mean((whitened-driftCoeff(1:numOfObservers(1),1))/whitened)*[1 1],[0 9],'r','LineWidth',2);
plot(mean((whitened-driftCoeff(1:numOfObservers(3),3))/whitened)*[1 1],[0 9],'b','LineWidth',2);

plot(mean((whitened-posCoeff(1:numOfObservers(2),2))/whitened)*[1 1],[0 9],'--k','LineWidth',2);
plot(mean((whitened-posCoeff(1:numOfObservers(1),1))/whitened)*[1 1],[0 9],'--r','LineWidth',2);
plot(mean((whitened-posCoeff(1:numOfObservers(3),3))/whitened)*[1 1],[0 9],'--b','LineWidth',2);

lh = legend(hh,{'drifts only','drifts + microsaccades'});
lh.FontSize = 10;
lh.EdgeColor = 'none';

xlabel('"Whitening" factor')
set(gca,'FontSize',14)
set(gca,'YTick',[3 5 7])
set(gca,'YTickLabel',{'Young','Age-Matched','AMD'})
xlim([0.2 1.2])


%% 
figure;
histogram(imCoeff,15);
xlabel('Power')
ylabel('Frequency')


%% eye movement spectra
%    'edriftF','epositionF','edPS','epPS','edPSse','epPSse'



figure;
subplot(1,2,1)
alp = 0.25;
th = 1.5;
PlotPatchedSE(edriftF(2:end),(edPS(2:end,2)),(edPSse(2:end,2)),'k',alp);
semilogx(edriftF(2:end),10*log10(edPS(2:end,2)),'-k','LineWidth',th)

hold on;
PlotPatchedSE(edriftF(2:end),(edPS(2:end,1)),(edPSse(2:end,1)),'r',alp);
semilogx(edriftF(2:end),10*log10(edPS(2:end,1)),'-r','LineWidth',th)
PlotPatchedSE(edriftF(2:end),(edPS(2:end,3)),(edPSse(2:end,3)),'b',alp);
semilogx(edriftF(2:end),10*log10(edPS(2:end,3)),'-b','LineWidth',th)
xlabel('Temporal frequency (hz)')
ylabel('Spectral density (dB)')
title('drifts only')
xlim([4 240])
grid on;
set(gca,'FontSize',14)

set(gca,'XTick',[10 100]);
set(gca,'XTickLabel',{'10','100'})

subplot(1,2,2)

PlotPatchedSE(edriftF(2:end),(epPS(2:end,2)),(epPSse(2:end,2)),'k',alp);
a(1) = semilogx(edriftF(2:end),10*log10(epPS(2:end,2)),'-k','LineWidth',th);
hold on;

PlotPatchedSE(edriftF(2:end),(epPS(2:end,1)),(epPSse(2:end,1)),'r',alp);
a(2) = semilogx(edriftF(2:end),10*log10(epPS(2:end,1)),'-r','LineWidth',th);

PlotPatchedSE(edriftF(2:end),(epPS(2:end,3)),(epPSse(2:end,3)),'b',alp);
a(3) = semilogx(edriftF(2:end),10*log10(epPS(2:end,3)),'-b','LineWidth',th);
xlabel('Temporal frequency (hz)')
ylabel('Spectral density (dB)')
lh= legend(a(1:3),{'Young','Age-Matched','AMD'});
lh.FontSize = 10;
xlim([4 240])
grid on;
title('drifts + miscrosaccades')
set(gca,'FontSize',14)
set(gca,'XTick',[10 100]);
set(gca,'XTickLabel',{'10','100'})



figure;
patientW = (whitened-driftCoeff(1:numOfObservers(3),3))/whitened;
% scatter(PRL, patientW); hold on;
tbl = table(PRL,patientW,'VariableNames',{'PRL','patientW'});
lm = fitlm(tbl,'patientW~PRL')
ph = lm.plot;
ph(1).MarkerSize = 10;
ph(1).Marker = 'd';
ph(1).MarkerEdgeColor = 'b';
ph(1).MarkerFaceColor = [.5 .5 1];
ph(2).Color = 'b';
ph(2).LineWidth = 2;
ph(3).Color = 'b';
ph(3).LineWidth = 2;
ph(4).Color = 'b';
ph(4).LineWidth = 2;
set(gca,'FontSize',14)
box off;
xlabel('PRL eccentricity (deg)')
ylabel('Whitening factor')
title('')

figure;
youngW = (whitened-driftCoeff(1:numOfObservers(2),2))/whitened;
oldW = (whitened-driftCoeff(1:numOfObservers(1),1))/whitened;
patientW = (whitened-driftCoeff(1:numOfObservers(3),3))/whitened;

fsAll = [fixStab(1:numOfObservers(2),2); fixStab(1:numOfObservers(1),1); fixStab(1:numOfObservers(3),3)];
wfAll = [youngW; oldW; patientW];

tbl2 = table(log10(fsAll),wfAll,'VariableNames',{'FixationStability','WhiteningFactor'});
lm2 = fitlm(tbl2,'WhiteningFactor~FixationStability')
scatter(log10(fixStab(1:numOfObservers(2),2)), youngW,msize,[.5 .5 .5],'filled','MarkerEdgeColor','k','Marker','o'); hold on;
scatter(log10(fixStab(1:numOfObservers(1),1)), oldW,msize,[1 .5 .5],'filled','MarkerEdgeColor','r','Marker','s');
scatter(log10(fixStab(1:numOfObservers(3),3)), patientW,msize,[.5 .5 1],'filled','MarkerEdgeColor','b','Marker','d');
plot(log10(linspace(min(fsAll),max(fsAll),100)),lm2.feval(log10(linspace(min(fsAll),max(fsAll),100))),'--k','LineWidth',2);
set(gca,'FontSize',14)
box off;
xlabel('log(Fixation stability) log(deg^2)')
ylabel('Whitening factor')
legend('Young','Age-Matched','AMD','Linear fit')
% xlim([-0.5 3.5])



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
