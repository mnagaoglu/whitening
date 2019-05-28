
clearvars
close all
clc;

load('diffusionData_n128.mat')

figure;
group = [diffusionData.group];
cols = [0 0 0; 1 0 0; 0 0 1];
ms = 100;
s = {'o','s','d'};

for i=[2 1 3]
    ix = group == i;
    scatter([diffusionData(ix).Dx_drift],[diffusionData(ix).Dy_drift],ms,...
        cols(i,:),s{i},'filled'); 
    hold on;
    
end

set(gca,'xscale','log','yscale','log','fontsize',14)
xlabel('D_x (arcmin^2/sec)')
ylabel('D_y (arcmin^2/sec)')
grid on
hold on; plot([.1 1000],[.1 1000],'--k')
xlim([1 1000])

legend('Young adults','Older adults','AMD patients')