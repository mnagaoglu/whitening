clearvars
% close all
clc

load(['../../ganglion/diffusionData_n128.mat']);
% load(['diffusionData_n128.mat']);

myColors = {'r','k','b'};
mySymbols = {'o','s','d'};
group = [diffusionData.group];

figure;
msize = 200;


for i=[2 1 3]
    ix = group == i;
    x = [diffusionData(ix).Dx_drift];
    y = [diffusionData(ix).Dy_drift];
    
    s = scatter(x,y,msize,myColors{i},mySymbols{i},'filled');  hold on;
    s.MarkerFaceAlpha = 0.7;
end

set(gca,'xscale','log','yscale','log','fontsize',14);
xlabel('D_x (arcmin^2/sec)')
ylabel('D_y (arcmin^2/sec)')
xlim([1 1000])
grid on;
plot([1 1000],[1 1000],'-','Color',[.7 .7 .7],'linewidth',2); hold on;
legend('Young adults','Older adults','AMD patients')
