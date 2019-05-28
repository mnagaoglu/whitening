% for FIGURE 1 in the paper
close all;
clearvars;
clc;

load('readyToAnalyzeData_new.mat');
groupNumber = [eyeMovements.group];

thePath = [pwd filesep 'Kumar and Chung 2014 reanalysis 2/'];
l = dir([thePath '*_sacsdrifts.mat']);

cols = [1 0 0; 0 0 0; 0 0 1];
figure;

toPlot = [25 52 62];
order = [2 1 3];
labels = {'Young','Age-matched','AMD'};
for i=1:length(l)
    if any(toPlot == i)
        posFile = [thePath num2str(i) '_480_hz_final_filtered.mat'];
        load(posFile,'eyePositionTraces','timeArray');

        
        subplot(1,3,order(toPlot == i))
        ix = timeArray < 10;
        plot(timeArray(ix), eyePositionTraces(ix,1),'-','linewidth',2,...
            'color',cols(groupNumber(i),:)); hold on;
        plot(timeArray(ix), eyePositionTraces(ix,2),':','linewidth',2,...
            'color',cols(groupNumber(i),:));
        set(gca,'fontsize',14);
%         if order(toPlot == i) == 3
            xlabel('Time (sec)')
%         end
%         if order(toPlot == i) == 2
            ylabel('Eye position (deg)')
%         end
        if order(toPlot == i) == 1
            legend('Hor','Ver','location','best');
        end
        ylim([-1. 1.5])
        xlim([0 5])
        grid off;
        box off;
        title(labels{order(toPlot == i)})

    end
    
end
