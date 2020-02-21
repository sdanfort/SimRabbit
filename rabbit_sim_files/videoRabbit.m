function videoRabbit(time,states,linkPosFun,fileName,targetPhi,targetQ,lc_,polyMod,xHat,plotFeatures)
frameRate = 30;
writerObj = VideoWriter(fileName,'MPEG-4');
writerObj.FrameRate = frameRate;
writerObj.Quality = 75;
open(writerObj);

plotPhiTarget = 0;
plotQTarget = 0;

if nargin > 4
    plotPhiTarget = 1;
end

if nargin > 5
    plotQTarget = 1;
end

if nargin < 7
    lc_ = [1,1,1];
end

if nargin < 8
    plotV = 0;
else
    plotV = 1;
end

if nargin < 10
    plotFeatures = [1,1,1];
end

ls_ = '-';

time = time-time(1);
% [time,i_] = unique(time);
% states = states(i_,:);

eventInds = find(~diff(time))+1;
eventInds = [1;eventInds(:);size(states,1)];

% [time,ind_] = unique(time);
% states = states(ind_,:);

link0 = linkPosFun(states(1,:).');
lstr = {};

h_fig = figure(20); clf

if plotV
    set(h_fig,'Position',[50,50,16*80,9*80])
    subplot(1,2,1);
    set(gca,'Position',[0.0427    0.11    0.4634    0.815])
else
    set(h_fig,'Position',[100,100,400,400])
end
hold on,box on
set(gca,'ytick',[])
xlabel('Position (m)')

if plotPhiTarget
    hPhi = plot(link0(3,1) + [0,-sin(targetPhi(1))*0.48],link0(3,2) + [0,cos(targetPhi(1))*0.48],'color',1 - 0.6*(1-lc_),'linestyle',':','linewidth',12);
    lstr = [lstr,'Pitch Target'];
end
if plotQTarget
    linkQ0 = linkPosFun(targetQ(1,:).');
    hq0 = plot(linkQ0(:,1),linkQ0(:,2),'color',1 - 0.2*(1-lc_),'linewidth',8,'linestyle','-'); %1 - 0.3*(1-[0,0.45,0.74])
    lstr = [lstr,'Configuration Target'];
end
hl = plot(link0(:,1),link0(:,2),'color',lc_,'linestyle',ls_,'linewidth',3);ha = gca;axis equal
lstr = [lstr,'Actual Configuration'];

plot(ha.XLim,[0,0],'k')

legend(lstr,'location','northwest')

yl = ha.YLim;
xl = ha.XLim;
drawnow

if plotV
    hAxV = subplot(1,2,2); hold on, box on
    set(hAxV,'position',[0.5652    0.1100    0.4109    0.8150])
    ssP = safeSetPlot(polyMod,80,80,2.5,0.05,30,hAxV,plotFeatures);
    ssP.update(xHat(1,3),xHat(1,4))
    hP = plot(xHat(1,1),xHat(1,2),'o','markersize',7,'linewidth',3,'markeredgecolor',0.4*lc_,'markerfacecolor',lc_);
    set(gca,'ylim',[-0.39,2.5])
    set(gca,'xlim',[polyMod.scaleParam.xMid{1}(1) - polyMod.scaleParam.xSize{1}(1),polyMod.scaleParam.xMid{8}(1) + polyMod.scaleParam.xSize{8}(1)])
end

tic
xOffset = link0(1,1);

currTime = 0;
frame(1,1) = getframe(h_fig);

for i = 1:length(eventInds)-1
    inds = eventInds(i):(eventInds(i+1)-1);
    while currTime < time(inds(end))
       st_ = interp1(time(inds),states(inds,:),currTime);
       link_ = linkPosFun(st_.');
       
       if plotPhiTarget
           targ_ = interp1(time(inds),targetPhi(inds),currTime);
           set(hPhi,'xdata',link_(3,1) + xOffset + [0,-sin(targ_)*0.48])
           set(hPhi,'ydata',link_(3,2) + [0,cos(targ_)*0.48])
       end
       if plotQTarget
           targQ0_ = interp1(time(inds),targetQ(inds,:),currTime);
           linkQ0_ = linkPosFun(targQ0_.');
           set(hq0,'xdata',linkQ0_(:,1)+xOffset);
           set(hq0,'ydata',linkQ0_(:,2));
       end
       if plotV
           xh_ = interp1(time(inds),xHat(inds,:),currTime);
           ssP.update(xh_(3),xh_(4))
           hP.XData = xh_(1);
           hP.YData = xh_(2);
           set(gca,'ylim',[-0.39,2.5])
           set(gca,'xlim',[polyMod.scaleParam.xMid{1}(1) - polyMod.scaleParam.xSize{1}(1),polyMod.scaleParam.xMid{8}(1) + polyMod.scaleParam.xSize{8}(1)])
       end
       
       set(hl,'xdata',link_(:,1)+xOffset);
       set(hl,'ydata',link_(:,2));
       set(ha,'ylim',yl);
       set(ha,'xlim',link_(3,1)+xOffset+xl)
       drawnow
       currTime = currTime + 1/frameRate;
       
       frame(1,end+1) = getframe(h_fig);
    end
    
    link_ = linkPosFun(states(inds(end),:).');
    xOffset = xOffset + link_(end,1);
end

writeVideo(writerObj,frame);
close(writerObj);

end
    