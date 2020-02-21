function videoSet(time,states,vFun,fileName,targetPhi,targetQ,lc_,ls_)
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

if nargin < 8
    ls_ = '-';
    lc_ = [1,1,1];
end

time = time-time(1);
% [time,i_] = unique(time);
% states = states(i_,:);

eventInds = find(~diff(time))+1;
eventInds = [1;eventInds(:);size(states,1)];

% [time,ind_] = unique(time);
% states = states(ind_,:);

link0 = linkPosFun(states(1,:).');
lstr = {};

h_fig = figure(20),clf,hold on

if plotPhiTarget
    hPhi = plot(link0(3,1) + [0,-sin(targetPhi(1))*0.48],link0(3,2) + [0,cos(targetPhi(1))*0.48],'color',1 - 0.6*(1-lc_),'linestyle',':','linewidth',6);
    lstr = [lstr,'Pitch Target'];
end
if plotQTarget
    linkQ0 = linkPosFun(targetQ(1,:).');
    hq0 = plot(linkQ0(:,1),linkQ0(:,2),'color',1 - 0.2*(1-lc_),'linewidth',4,'linestyle','-'); %1 - 0.3*(1-[0,0.45,0.74])
    lstr = [lstr,'Configuration Target'];
end
hl = plot(link0(:,1),link0(:,2),'color',lc_,'linestyle',ls_,'linewidth',1.5);ha = gca;axis equal
lstr = [lstr,'Actual Configuration'];

hold on,plot(ha.XLim,[0,0],'k')

legend(lstr,'location','northwest')

yl = ha.YLim;
xl = ha.XLim;
drawnow

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
    