function rabbitStillsTime(time,states,linkPosFun,nStills,tScale,targetPhi,targetQ,lc_)

plotPhiTarget = 0;
plotQTarget = 0;
if nargin > 5
    plotPhiTarget = 1;
end

if nargin > 6
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
eventInds = [1;eventInds(:)];
time(eventInds) = time(eventInds) + 1e-9;

link0 = linkPosFun(states(1,:).');
lstr = {};
figure,clf,hold on
if plotPhiTarget
    plot(link0(3,1) + [0,-sin(targetPhi(1))*0.48],link0(3,2) + [0,cos(targetPhi(1))*0.48],'color',1 - 0.4*(1-lc_),'linestyle',':','linewidth',5);
    lstr = [lstr,'Pitch Target'];
end
if plotQTarget
    linkQ0 = linkPosFun(targetQ(1,:).');
    plot(linkQ0(:,1),linkQ0(:,2),'color',1 - 0.4*(1-lc_),'linewidth',3,'linestyle','-'); %1 - 0.3*(1-[0,0.45,0.74])
    lstr = [lstr,'Configuration Target'];
end
plot(link0(:,1) - link0(3,1),link0(:,2),'color',lc_,'linestyle','-','linewidth',1.5);
lstr = [lstr,'Actual Configuration'];

ha = gca;axis equal, hold on,


tVals = cell(nStills+1,1);
tVals{1} = '0';

for i = 1:nStills
    currTime = i/nStills*time(end);
    tVals{i+1} = sprintf('%3.1f', currTime);
        
    st_ = interp1(time,states,currTime);
    link_ = linkPosFun(st_.');
    
    xHip = link_(3,1);
    
    if plotPhiTarget
        targ_ = interp1(time,targetPhi,currTime);
        plot(currTime*tScale + [0,-sin(targ_)*0.48],link_(3,2) + [0,cos(targ_)*0.48],'color',1 - 0.4*(1-lc_),'linestyle',':','linewidth',5);
    end
    if plotQTarget
        targQ0 = interp1(time,targetQ,currTime);
        linkQ0 = linkPosFun(targQ0.');
        plot(linkQ0(:,1)  - xHip + currTime*tScale,linkQ0(:,2),'color',1 - 0.4*(1-lc_),'linewidth',3,'linestyle','-'); %1 - 0.3*(1-[0,0.45,0.74])
    end
    plot(link_(:,1) - xHip + currTime*tScale ,link_(:,2),'color',lc_,'linestyle','-','linewidth',1.5);
end

plot(ha.XLim,[0,0],'k','linewidth',2)
set(gca,'xtick',linspace(0,tScale*time(end),nStills+1),'xticklabels',tVals)
xlabel('Time (s)')
% legend(lstr,'location','northwest')

end