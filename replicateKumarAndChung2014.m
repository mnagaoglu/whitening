%% save data in a format that ReVAS can chew on.
% close all;
% clearvars;
% clc;
% 
% thePath = '/Users/mnagaoglu/Desktop/data/';
% load('readyToAnalyzeData_new.mat');
% 
% try
%     for i=1:length(eyeMovements)
%         %if eyeMovements(i).group == 1
%             timeArray = eyeMovements(i).rawTime';
%             eyePositionTraces = eyeMovements(i).rawPosition;
%             parametersStructure = struct;
%             referenceFramePath = '';
%             subject = eyeMovements(i).initials;
% 
%             save([thePath num2str(i) '_480_hz_final.mat'],...
%                 'eyePositionTraces','timeArray','referenceFramePath',...
%                 'parametersStructure','subject');
%         %end
%     end
% 
% catch err
%     err;
% end
% 

%% read filtered data and parsed saccades and plot them
% 
% close all;
% clearvars;
% clc;
% 
% thePath = '/Users/mnagaoglu/Desktop/data/';
% l = dir([thePath '*_sacsdrifts.mat']);
% 
% figure;
% for i=24:length(l)
%     saccadeFile = [thePath num2str(i) '_480_hz_final_filtered_sacsdrifts.mat'];
%     posFile = [thePath num2str(i) '_480_hz_final_filtered.mat'];
%     load(saccadeFile,'saccades');
%     load(posFile,'eyePositionTraces','timeArray');
%     
%     cla;
%     plot(timeArray,eyePositionTraces(:,1),'-r','linewidth',2); hold on;
%     plot(timeArray,eyePositionTraces(:,2),'-b','linewidth',2); 
%     for j=1:length(saccades)
%         ix = saccades(j).onsetIndex:saccades(j).offsetIndex;
%         plot(timeArray(ix),eyePositionTraces(ix,1),'or','markersize',10);
%         plot(timeArray(ix),eyePositionTraces(ix,2),'ob','markersize',10);
%     end
%     legend('hor','ver');
%     xlabel('time (sec)')
%     ylabel('eye position (deg)')
%     
% end



%% count microsaccades and compute avg size

close all;
clearvars;
clc;

thePath = '/Users/mnagaoglu/Desktop/data/';
l = dir([thePath '*_sacsdrifts.mat']);
load('readyToAnalyzeData_new.mat');


% read parsed saccades and compile a list of their statistics
amplitudes = [];
counts = [];
maxTime = [];
% initials = [];

for i=1:length(l)
    saccadeFile = [thePath num2str(i) '_480_hz_final_filtered_sacsdrifts.mat'];
    posFile = [thePath num2str(i) '_480_hz_final_filtered.mat'];
    load(saccadeFile,'saccades');
    load(posFile,'eyePositionTraces','timeArray');
    
%     originalFile = [thePath num2str(i) '_480_hz_final.mat'];
%     load(originalFile,'subject');
%     initials = [initials; subject];
    
    % remove saccades larger than 5deg
    saccades = saccades([saccades.vectorAmplitude] < 3.0);
    
    maxTime = [maxTime; 10 max(timeArray)];
    firstTenSeconds = [saccades.onsetTime] < 10;
    counts = [counts; sum(firstTenSeconds) length(saccades)];
    amplitudes = [amplitudes; ...
        mean([saccades(firstTenSeconds).vectorAmplitude])...
        mean([saccades.vectorAmplitude])]; %#ok<*AGROW>
    
end

% take within-subject averages
initials = lower({eyeMovements.initials});
% [subjects, ia, ic] = unique(initials);
groupNumber = [eyeMovements.group];

rate = [];
avgAmp = [];
group = [];
subjects = [];
for i=1:3
    ix = find(i == groupNumber); 
    
    thisGroupInitials = initials(ix);
    [thisGroupSubjects, ia, ic] = unique(thisGroupInitials);
    subjects = [subjects; thisGroupSubjects'];
    for j=1:length(thisGroupSubjects)
        indices = ix(j == ic);
        rate = [rate; mean(counts(indices,:)./maxTime(indices,:),1)];
        avgAmp = [avgAmp; mean(amplitudes(indices,:),1)];
        group = [group; i];
    end

end


%% plot results
figure('units','normalized','outerposition',[.5 .5 .4 .3]);
d = 0.2;
cols = [0 0 0; 1 0 0; 0 0 1];
ms = 10;
order = [2 1 3];
for i=1:3
    ix = find(group == order(i));
    
    subplot(1,2,1)
    p1 = plot(i-d,rate(ix,1),'o-','Color',cols(i,:),...
        'markersize',ms); hold on;
    p2 = plot(i+d,rate(ix,2),'o-','Color',cols(i,:),...
        'markerfacecolor',cols(i,:),'markersize',ms);
    for j=1:length(ix)
        plot(i+[-d d],rate(ix(j),:),'-','Color',cols(i,:),'linewidth',2);
%         text(i-d-0.35, rate(ix(j),1), subjects{ix(j)},'color',cols(i,:),...
%             'fontsize',14);
    end
    xlim([0 4]);
    set(gca,'fontsize',14,'xtick',[1 2 3],'xticklabel',{'Young','Old','AMD'})
    if i ==1
        %legend([p1; p2],{'<10sec';'all'});
        title('Microsaccade rate')
    end
    
    subplot(1,2,2)
    plot(i-d,avgAmp(ix,1),'o-','Color',cols(i,:),...
        'markersize',ms); hold on;
    plot(i+d,avgAmp(ix,2),'o-','Color',cols(i,:),...
        'markerfacecolor',cols(i,:),'markersize',ms);
    for j=1:length(ix)
        plot(i+[-d d],avgAmp(ix(j),:),'-','Color',cols(i,:),'linewidth',2);
%         text(i-d-0.35, avgAmp(ix(j),1), subjects{ix(j)},'color',cols(i,:),...
%             'fontsize',14);
    end
    xlim([0 4]);
    set(gca,'fontsize',14,'xtick',[1 2 3],'xticklabel',{'Young','Old','AMD'})
    if i ==1
        title('Average amplitude (deg)')
    end
    
end




%% compute and plot within-group averages
figure('units','normalized','outerposition',[.5 .5 .3 .3]);
d = 0.2;
cols = [0 0 0; 1 0 0; 0 0 1];
ms = 10;
order = [2 1 3];

groupAvgRate = [];
groupAvgAmp = [];
groupSERate = [];
groupSEAmp = [];
for i=1:3
    ix = find(group == order(i));
    groupAvgRate(i,:) = mean(rate(ix,:),1);
    groupAvgAmp(i,:) = mean(avgAmp(ix,:),1);
    groupSERate(i,:) = groupAvgRate(i,:)/sqrt(length(ix));
    groupSEAmp(i,:) = groupAvgAmp(i,:)/sqrt(length(ix));
    
    subplot(1,2,1)
    p1 = plot(i-d,groupAvgRate(i,1),'o-','Color',cols(i,:),...
        'markersize',ms); hold on;
    p2 = plot(i+d,groupAvgRate(i,2),'o-','Color',cols(i,:),...
        'markerfacecolor',cols(i,:),'markersize',ms);
    errorbar(i+[-d d],groupAvgRate(i,:),groupSERate(i,:),'-',...
        'Color',cols(i,:),'linewidth',2);

    xlim([0 4]);
    set(gca,'fontsize',14,'xtick',[1 2 3],'xticklabel',{'Young','Old','AMD'})
    if i ==1
        %legend([p1; p2],{'<10sec';'all'});
        title('Microsaccade rate')
    end
    
    subplot(1,2,2)
    p1 = plot(i-d,groupAvgAmp(i,1),'o-','Color',cols(i,:),...
        'markersize',ms); hold on;
    p2 = plot(i+d,groupAvgAmp(i,2),'o-','Color',cols(i,:),...
        'markerfacecolor',cols(i,:),'markersize',ms);
    errorbar(i+[-d d],groupAvgAmp(i,:),groupSEAmp(i,:),'-',...
        'Color',cols(i,:),'linewidth',2);
    
    xlim([0 4]);
    set(gca,'fontsize',14,'xtick',[1 2 3],'xticklabel',{'Young','Old','AMD'})
    if i ==1
        %legend([p1; p2],{'<10sec';'all'});
        title('Average amplitude (deg)')
    end
    
end

