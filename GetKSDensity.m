function [density, X, Y, bandwidth, PRL, fh, stats] = GetKSDensity(xDeg, yDeg, isShowPlot)

if nargin<3
    isShowPlot = 0;
end
data = [xDeg yDeg];

use = sum(abs(data),2) < 100;
xy = [xDeg(use) yDeg(use)];

% Perform 2d kernel density estimation on the gaze data
kdetrim = 0.045; % set a cutoff for the trimming function
bw = std(xy)/length(xy)^(1/6); % bandwidth size
[bandwidth,density,X,Y] = kde2d_trimmed(kdetrim, xy, 128,[], [], bw);

if isnan(density)
    density = NaN;
    X = NaN;
    Y = X;
    bandwidth = [];
    PRL = [];
    fh = [];
    stats = [];
    fprintf('Not sufficient data for kernel density estimation!!!!\n');
    return;
end

maxdensity = max(density(:));
I = find(maxdensity == density);
PRL = [X(I), Y(I)];




fh = figure;
h = polar([0 2*pi], [0 15]);
delete(h);
hold on;

msize = 20;
mcolor = 'k';
sh = scatter(xDeg, yDeg,msize,mcolor,'filled','MarkerEdgeColor','none');
sh.MarkerFaceAlpha = 0.3;

steps = linspace(min(density(:)),max(density(:)),19);
steps = steps(2:end-1);
[C, ch] = contour(X,Y,reshape(density,size(X,1),size(X,2)),steps); hold on
plot(PRL(1),PRL(2),'+r','MarkerSize',15);
caxis([0 .075]);
colormap(jet)




% STATS: compute mean and dispersion (sqrt(mean of variacne)), BCEA and overlay the principal axes
bceathd = 0.68; % 68% CI
mu = bimean(X,Y,density);
[pv, pd] = bivar(X,Y,density);
dispersion = sqrt(mean(pv));
bcea = pi*chi2inv(bceathd,2)*sqrt(prod(pv)); % BCEA



stats.mu = mu;
stats.pv = pv;
stats.pd = pd;
stats.dispersion = dispersion;
stats.bcea = bcea;



%% isoline method
cumProb = sum(density(:));
isoAreas = [];
isoProb = [];
for i=1:length(steps)

    try
        if i==1
            st = 2;
            en = C(2,1)+1;
        else
            st = en+2;
            en = C(2,st-1)+st-1;
        end
        
        islands = find(C(1,:) == steps(i));
        in = false(length(X(:)),1);
        tempArea = 0;
        for j=1:length(islands)   
            st = islands(j)+1;
            xv = C(1,st:C(2,st-1)+st-1);
            yv = C(2,st:C(2,st-1)+st-1);
            in = in | inpolygon(X(:),Y(:),xv',yv');
            tempArea = tempArea + polyarea(xv,yv);
        end

        isoProb(i) = sum(density(in))/cumProb;
        isoAreas(i) = tempArea;
        
    catch err
        err.message
        err.stack.line
        err.stack.name
        fprintf('Error occured during isoline method\n');
        rethrow(err)
    end
end


try
    desiredStep = interp1(isoProb,steps,0.68,'pchip');
    figure(fh);
    [C, ch] = contour(X,Y,reshape(density,size(X,1),size(X,2)),[desiredStep 1]);
    delete(ch);


    islands = find(C(1,:) == desiredStep);
    isoa = 0;
    for i=1:length(islands)
        st = islands(i)+1;
        xv = C(1,st:C(2,st-1)+st-1);
        yv = C(2,st:C(2,st-1)+st-1);
        isoa = isoa + polyarea(xv,yv);
    end
    stats.areaInContour = isoa;
    stats.isoProb = isoProb;
    stats.isoAreas = isoAreas;
catch erriso
    erriso.message
    fprintf('\n\nError during interpolation or isoline area computation...\n')
    rethrow(erriso)
end


if isShowPlot
    figure(fh);
    title(sprintf('PRL %.2f %.2f \n BCEA: %.2f, ISOA: %.2f',PRL(1), PRL(2),...
        bcea, isoa));
    set(gca,'fontsize',16)
else
    if ~isempty(fh)
        delete(fh);
    end
    fh = [];
end


