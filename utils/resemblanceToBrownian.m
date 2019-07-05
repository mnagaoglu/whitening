clearvars
close all
clc

load('resemblance.mat')

load('diffusionData_n128.mat')
n = 3;

myColors = {'r','k','b'};
mySymbols = {'o','s','d'};
group = [diffusionData.group];
initials = {diffusionData.initials};
msize = 100;
Fs = 480;
N = 128;
t = 0:1/Fs:(N-1)/Fs; 
alp = 0.5;

figure; hold on;
yloc = [2 1 3];
for i=1:3
    ix = find(group == i);
    gofs = [diffusionData(ix).gof_drift];
    
    
%     % compute rmse
%     rhat = nan(size(r));
%     for j=1:length(ix)
%         temp = diffusionData(ix(j)).D_fit_drift * 4 * t(2:end);
%         rhat(j) = rms(temp - diffusionData(ix(j)).meanD2_drift);
%     end
    
    r = [gofs.rsquare];
    c = correlation(n,ix,3);
    rmse = rootMeanSquared(n,ix,3);
    rnd = randn(length(ix),1);
    
    subplot(3,1,1)
    s(i) = scatter(r,yloc(i) + 0.025*rnd,msize,myColors{i},mySymbols{i},'filled'); 
    s(i).MarkerFaceAlpha = alp;
    hold on;
    grid on;
    set(gca,'fontsize',14,'xscale','linear') 
    ylim([0 4]);
    xlim([0.75 1.01]);
    set(gca,'ytick',[1 2 3],'yticklabel',{'young','older','MD'})
    xlabel('R^2')
    
    subplot(3,1,2)
    s2(i) = scatter(rmse,yloc(i) + 0.025*rnd,msize,myColors{i},mySymbols{i},'filled'); 
    s2(i).MarkerFaceAlpha = alp;
    hold on;
    set(gca,'fontsize',14,'xscale','log')
    ylim([0 4]);
    grid on;
    xlim([0.01 100]);
    set(gca,'ytick',[1 2 3],'yticklabel',{'young','older','MD'})
    xlabel('RMS error')
    
    subplot(3,1,3)
    s3(i) = scatter(c,yloc(i) + 0.025*rnd,msize,myColors{i},mySymbols{i},'filled'); 
    s3(i).MarkerFaceAlpha = alp;
    hold on;
    grid on;
    set(gca,'fontsize',14,'xscale','linear')
    ylim([0 4]);
    xlim([0.96 1.01]);
    set(gca,'ytick',[1 2 3],'yticklabel',{'young','older','MD'})
    xlabel('Pearson''s \rho')
end

% SaveAsPDF(gcf,'gof.pdf')



figure('units','normalized','outerposition',[.3327    0.6381    0.4863    0.3333]); 
hold on;
plot([.001 1000],[.001 1000],'-','color',[1 1 1]*.7,'linewidth',2);
for i=[2 1 3]
    ix = find(group == i);
    gofs = [diffusionData(ix).gof_drift];
    D = [diffusionData(ix).D_fit_drift];
    r = [gofs.rsquare];
    c = correlation(n,ix,3);
    rmse = rootMeanSquared(n,ix,3);
    rnd = randn(length(ix),1);
    
%     % compute rmse
%     rhat = nan(size(r));
%     for j=1:length(ix)
%         temp = diffusionData(ix(j)).D_fit_drift * 4 * t(2:end);
%         rhat(j) = rms(temp - diffusionData(ix(j)).meanD2_drift);
%     end

  
    subplot(1,3,1)
    s(i) = scatter(r,D,msize,myColors{i},mySymbols{i},'filled'); 
    s(i).MarkerFaceAlpha = alp;
    hold on;
    grid on;
    set(gca,'fontsize',14,'xscale','linear','yscale','log') 
    xlim([0.75 1.01]);
    xlabel('R^2')
    ylim([1 1000])
    set(gca,'ytick',[1 10 100 1000],'yticklabel',[1 10 100 1000])
    ylabel('D (arcmin^2/sec)')
    
    subplot(1,3,2)
    s2(i) = scatter(rmse,D,msize,myColors{i},mySymbols{i},'filled'); 
    s2(i).MarkerFaceAlpha = alp;
    hold on;
    set(gca,'fontsize',14,'xscale','log','yscale','log',...
        'xtick',[.01 .1 1 10 100],'xticklabel',[.01 .1 1 10 100])
    grid on;
    xlim([0.01 100]);
    set(gca,'ytick',[1 10 100 1000],'yticklabel',[1 10 100 1000])
    xlabel('RMS error')
    ylim([1 1000])
    
    subplot(1,3,3)
    s3(i) = scatter(c,D,msize,myColors{i},mySymbols{i},'filled'); 
    s3(i).MarkerFaceAlpha = alp;
    hold on;
    grid on;
    set(gca,'fontsize',14,'xscale','linear','yscale','log')
    set(gca,'ytick',[1 10 100 1000],'yticklabel',[1 10 100 1000])
    xlim([0.96 1.001]);
    xlabel('Pearson''s \rho')
    ylim([1 1000])
end



subplot(1,3,2)
legend({'young','older','MD'},'location','best')

% SaveAsPDF(gcf,'gof_scatter.pdf')
