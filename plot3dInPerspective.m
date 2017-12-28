function plot3dInPerspective(sf,tf,S, ...
    isColorbar,isX,isY,titlestr,cols)

if nargin<4
    isColorbar = 0;
end

newSf = logspace(min(log10(sf)),max(log10(sf)),512);
newTf = logspace(min(log10(tf)),max(log10(tf)),512);
[X,Y] = meshgrid(sf,tf);
[Xq,Yq] = meshgrid(newSf,newTf);
newS = interp2(X,Y,S,Xq,Yq,'cubic');


surf(newSf,newTf,newS,'EdgeColor','none','LineWidth',0.25,'EdgeAlpha',.5);
xlim([.25 65])
ylim([4 240])
minZ = -5;
maxZ = 72;
zlim([minZ maxZ])
view(2);set(gca,'XScale','log');
set(gca,'YScale','log');
colormap(jet(256))
set(gca,'CLim',[minZ maxZ])
set(gca,'XTick',[1 3 10 30])
set(gca,'XTickLabel',{'1','3','10','30'})
set(gca,'YTick',[3 10 30 100]);
set(gca,'YTickLabel',{'3','10','30','100'})
ax = gca;
ax.TickDir = 'out';
ax.TickLength = [0.02 0.02];
title(titlestr);

if isX
    xlabel('Spatial frequency (cpd)');
else
    xlabel(' ')
end
if isY
	ylabel('Temporal frequency (Hz)');
else
    ylabel(' ')
end
set(gca,'FontSize',12);
pause(.2);
s1Pos = get(gca,'position');
if isColorbar
    cb = colorbar('eastoutside');
    cb.Label.String = 'Spectral Density (dB)';
    cb.Color = cols;
end
pause(.2);
s2Pos = get(gca,'position');
s2Pos(3:4) = s1Pos(3:4);
set(gca,'position',s2Pos);

if isColorbar
   cb.Position = cb.Position./[1 1 3 1]; 
end


