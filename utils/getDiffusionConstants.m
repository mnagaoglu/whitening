function getDiffusionConstants(N)
%% getDiffusionConstants
% This function loads the pre-processed, ready-to-analyze eye movement
% traces and estimates diffusion coefficients by assuming that FEM follow a
% standard isotropic or anisotropic Brownian motion.
%
% We compute diffusion constants for different temporal window lengths..
% just to see...

try
    load('organizedData.mat','organizedData');
catch err
    err.message
    err.stack.line
    err.name
    disp('organizedData.mat cannot be found.')
    return;
end

if nargin<1
    N = [32 64 128 256 512 1024 2048];
end

rootMeanSquared = nan(length(N),length(organizedData),3);
correlation = rootMeanSquared;
rsq = correlation;

for j=1:length(N)
    
    n = N(j); % number of samples (window width: Fs*n seconds) 

    % compute diffusion constant for each trial
    for i=1:length(organizedData)
        ct = organizedData(i);

        % drifts only
        [newCt.D_drift, newCt.Dx_drift, newCt.Dy_drift, newCt.t_drift,...
            newCt.D_se_drift, newCt.Dx_se_drift, newCt.Dy_se_drift,...
            newCt.meanD2x_drift, newCt.seD2x_drift, newCt.meanD2y_drift,...
            newCt.seD2y_drift, newCt.meanD2_drift, newCt.seD2_drift,...
            newCt.Dx_fit_drift, newCt.Dy_fit_drift,newCt.D_fit_drift, newCt.CIx_drift,...
            newCt.CIy_drift, newCt.CI_drift, ...
            newCt.gofx_drift, newCt.gofy_drift, newCt.gof_drift] = ...
            estimateDiffusionConstant(ct.stitchedPosition,480,n);

        % drifts + microsaccades
        [newCt.D, newCt.Dx, newCt.Dy, newCt.t, newCt.D_se, newCt.Dx_se, newCt.Dy_se,...
            newCt.meanD2x, newCt.seD2x, newCt.meanD2y, newCt.seD2y, newCt.meanD2, newCt.seD2,...
            newCt.Dx_fit, newCt.Dy_fit,newCt.D_fit, newCt.CIx, newCt.CIy, newCt.CI,...
            newCt.gofx, newCt.gofy, newCt.gof] = ...
            estimateDiffusionConstant(ct.newPos,480,n);

        newCt.initials = ct.initials;
        newCt.fullFileName = ct.fullFileName;
        newCt.group = ct.group;

        diffusionData(i) = newCt;
        
        x = newCt.Dx_fit_drift*2*newCt.t_drift/3600;
        y = newCt.Dy_fit_drift*2*newCt.t_drift/3600;
        v = newCt.D_fit_drift*4*newCt.t_drift/3600;

%         figure(1);
%         cla;
%         p1 = plot(newCt.t_drift,newCt.meanD2x_drift*3600,'-','LineWidth',2); hold on;
%         p2 = plot(newCt.t_drift,newCt.meanD2y_drift*3600,'-','LineWidth',2);
%         p3 = plot(newCt.t_drift,newCt.meanD2_drift*3600,'-','LineWidth',2);
%         plot(newCt.t_drift,x*3600,'--',...
%             'LineWidth',2,'color',p1.Color); hold on;
%         plot(newCt.t_drift,y*3600,'--',...
%             'LineWidth',2,'color',p2.Color);
%         plot(newCt.t_drift,v*3600,'--',...
%             'LineWidth',2,'color',p3.Color);  
% %         plot(newCt.t_drift,x,':',...
% %             'LineWidth',2,'color',p1.Color); hold on;
% %         plot(newCt.t_drift,y,':',...
% %             'LineWidth',2,'color',p2.Color);
% %         plot(newCt.t_drift,v,':',...
% %             'LineWidth',2,'color',p3.Color);
%         set(gca,'fontsize',14) 
%         xlabel('Time lag (sec)')
%         ylabel('\Deltap^2 (arcmin^2)')
%         title(num2str(newCt.gof_drift.rmse))
        
        
        % look at the RMS error as resemblance
        rootMeanSquared(j,i,1) = 3600*rms(newCt.meanD2x_drift-x);
        rootMeanSquared(j,i,2) = 3600*rms(newCt.meanD2y_drift-y);
        rootMeanSquared(j,i,3) = 3600*rms(newCt.meanD2_drift-v);
        
        correlation(j,i,1) = corr(newCt.meanD2x_drift',x');
        correlation(j,i,2) = corr(newCt.meanD2y_drift',y');
        correlation(j,i,3) = corr(newCt.meanD2_drift',v');
        
        rsq(j,i,1) = newCt.gofx_drift.rsquare;
        rsq(j,i,2) = newCt.gofy_drift.rsquare;
        rsq(j,i,3) = newCt.gof_drift.rsquare;
        

    end

    save(sprintf('diffusionData_n%d.mat',n),'diffusionData');

end

save('resemblance.mat','rootMeanSquared','correlation','rsq');



function [D, Dx, Dy, t, D_se, Dx_se, Dy_se, ...
    meanD2x, seD2x, meanD2y, seD2y, meanD2, seD2,...
    Dx_fit, Dy_fit,D_fit, CIx, CIy, CI, gofx, gofy, gof] = ...
    estimateDiffusionConstant(pos,Fs,n)

overlap = 0.5; % fraction of overlap between consecutive windows
st = 1:round(n*(1-overlap)):(length(pos)-n); % start indices of each window
t = 0:1/Fs:(n-1)/Fs; % time axis for plotting

timeLimit = 30; % seconds
maxNumberOfSamples = round(timeLimit*Fs);
if maxNumberOfSamples < length(pos)
%     pos = pos(1:maxNumberOfSamples,:);
end

try
    for i=1:length(st)
        % subtract the initial point so that all traces start at 0
        Wx(i,:) = ((pos(st(i):st(i)+n-1,1)) - pos(st(i),1))'; 
        Wy(i,:) = ((pos(st(i):st(i)+n-1,2)) - pos(st(i),2))';
    end
catch err
    err.message
    err.stack.line
    err.stack.name
end
% Compute displacement-squared
[meanD2x, seD2x, meanD2y, seD2y, meanD2, seD2] = ...
    DisplacementSquared(Wx,Wy,length(st));

% Estimate Diffusion constants for anisotropic Brownian motion
% Note the additional factor of 2 for two-dimensional case to account for
% larger number of dimensions.
t = t(2:end);
Dx = nanmean(meanD2x./(2*t))*3600;
Dy = nanmean(meanD2y./(2*t))*3600;
D = nanmean(meanD2./(2*2*t))*3600; 

Dx_se = nanmean(seD2x./(2*t))*3600;
Dy_se = nanmean(seD2y./(2*t))*3600;
D_se = nanmean(seD2./(2*2*t))*3600; 

% estimate diffusion by fitting a line, where each displacement squared is
% downweighted with its standard error
[Dx_fit,~, CIx, gofx] = fitTimeVsDisplacementSquared(t, meanD2x*3600, 1./seD2x, 1);
[Dy_fit,~, CIy, gofy] = fitTimeVsDisplacementSquared(t, meanD2y*3600, 1./seD2y, 1);
[D_fit, fitresult, CI, gof]  = fitTimeVsDisplacementSquared(t, meanD2*3600, 1./seD2, 2);


% figure(12312); cla;
% scatter(t, 3600*meanD2,150,'filled'); hold on;
% plot(t,feval(fitresult,t),'-','linewidth',2);
% title([num2str(gof.rmse) ' ' num2str(gof.rmse2)])
% set(gca,'fontsize',12);


function [d, fitresult, ci, gof] = fitTimeVsDisplacementSquared(x, y, w, ndim)

% prepare data
[xData, yData, weights] = prepareCurveData( x, y, w ); 

% % normalize weights to sum to 1
weights = weights / sum(weights);

% weights = ones(size(weights));

% Set up fittype and options.
ft = fittype( 'a*x', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = 0.45733407813134;
opts.Weights = weights;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% ci = predint(fitresult,x,0.95,'functional','on');
ci = confint(fitresult)/(2*ndim);
d = fitresult.a/(2*ndim);

gof.rmse2 = rms(feval(fitresult,xData) - yData);




% A helper function for computing mean and standard error of 
% displacement-squared, i.e., $\Delta D^{2}$ 
function [meanD2x, seD2x, meanD2y, seD2y, meanD2, seD2] = ...
    DisplacementSquared(Wx,Wy,nsamples)

% Displacement-squared 
D2x = Wx(:,2:end).^2;
D2y = Wy(:,2:end).^2;
D2 = Wx(:,2:end).^2 + Wy(:,2:end).^2;

% means
meanD2x = nanmean(D2x,1);
meanD2y = nanmean(D2y,1);
meanD2 = nanmean(D2,1);

% standard error
seD2x = nanstd(D2x,[],1)/sqrt(nsamples);
seD2y = nanstd(D2y,[],1)/sqrt(nsamples);
seD2 = nanstd(D2,[],1)/sqrt(nsamples);
