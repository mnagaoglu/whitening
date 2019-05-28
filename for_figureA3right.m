clearvars
close all
clc

load(['../ganglion/diffusionData_n128.mat']);


myColors = {'r','k','b'};
mySymbols = {'o','s','d'};
group = [diffusionData.group];

figure;
msize = 200;
plot([1 1000],[1 1000],'-','Color',[.7 .7 .7],'linewidth',2); hold on;
for i=1:3
    ix = group == i;
    x = [diffusionData(ix).Dx_drift];
    y = [diffusionData(ix).Dy_drift];
    
    s(i) = scatter(x,y,msize,myColors{i},mySymbols{i},'filled'); 
    s(i).MarkerFaceAlpha = 0.7;
end

set(gca,'xscale','log','yscale','log','fontsize',20);
xlabel('D_x (arcmin^2/sec)')
ylabel('D_y (arcmin^2/sec)')
xlim([1 1000])

legend(s,{'Young adults','Older adults','AMD patients'})
