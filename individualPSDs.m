clearvars
close all
clc


load('FinalResultsClean.mat');
load('WF.mat');
wf = {youngWF; agedWF; patientWF};
% wf = {youngWF_pos; agedWF_pos; patientWF_pos};

cols = [0 0 0; 1 0 0 ; 0 0 1];
labels = {'Young adults','Older adults','AMD patients'};
th = 1.5;
figure('units','normalized','outerposition',[0.2000    0.2600    0.5286    0.3400]);
groupNo = [2 1 3]; % young, old, amd
band = cell(3,1);
for i=1:3
    group = groupNo(i);
    subplot(1,3,i)
    thisPSD = squeeze(nanmean(10*log10(drift_psd(2:end,2:end,:,group)),1));
    semilogx(sf(2:end),thisPSD,'-','Color',cols(i,:),'LineWidth',th); hold on;
    grid on;
    set(gca,'fontsize',14,'xtick',[1 10]);
    if i == 2
        xlabel('Spatial frequency (cpd)')
    end
    
    if i==1
        ylabel('Power (dB)');
    end
    xlim([.2 64])
    ylim([0 40])
    title(labels{i});
    
    % get the 3dB cutoff of the whitened region
    try
    for j=1:size(thisPSD,2)
        cutoff = sf(find(thisPSD(:,j) < (max(thisPSD(:,j)) - 3),1,'first'));
        if ~isempty(cutoff)
            band{i} = [band{i}; cutoff];
        end
    end
    catch err
        err;
    end
    
end



% plot the scattergram of whitening factors and whitened bandwidths
figure;
symbols = {'o','s','d'};
ms = 100;
for i=1:3
    scatter(band{i}, wf{i},ms,cols(i,:),symbols{i},'filled'); hold on;
end
set(gca,'fontsize',14,'xscale','log','yscale','linear');
box off;
xlabel('Whitening 3dB cutoff (cpd)')
ylabel('Whitening factor');
legend('Yound adults','Older adults','AMD patients','location','best');
xlim([.5 50]);
ylim([0 1.1])
set(gca,'xtick',[1 10])
grid on

