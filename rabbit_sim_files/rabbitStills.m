function rabbitStills(time,states,linkPosFun,nStills,xSpace,targetPhi,targetQ,lc_)
if nargin < 5
    addSpace = 0;
else
    addSpace = 1;
end
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
time(eventInds) = time(eventInds) + 1e-6;
% eventInds = [1;eventInds(:);size(states,1)];

% [time,ind_] = unique(time);
% states = states(ind_,:);

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
plot(link0(:,1),link0(:,2),'color',lc_,'linestyle','-','linewidth',1.5);
lstr = [lstr,'Actual Configuration'];

tVals = cell(nStills+1,1);
hipX = zeros(nStills+1,1);
tVals{1} = 0;
hipX(1) = link0(3,1);

ha = gca;axis equal, hold on,

xOffset = link0(1,1);
for i = 2:length(eventInds)
    link_ = linkPosFun(states(eventInds(i)-1,:).');
    xOffset(i) = xOffset(i-1) + link_(end,1);
end

xOff1_ = 0;
for i = 1:nStills
    currTime = i/nStills*time(end);
    tVals{i+1} = sprintf('%3.1f', currTime);

    evInd_ = find(currTime+1e-6 <= [time(eventInds);inf],1);
    xOff2_ = xOffset(evInd_-1);
    
    if addSpace
        xOff1_ = xOff1_ + xSpace;
    end
    
    st_ = interp1(time,states,currTime);
       
    link_ = linkPosFun(st_.');
    if plotPhiTarget
        targ_ = interp1(time,targetPhi,currTime);
        plot(link_(3,1) + xOff1_ + xOff2_ + [0,-sin(targ_)*0.48],link_(3,2) + [0,cos(targ_)*0.48],'color',1 - 0.4*(1-lc_),'linestyle',':','linewidth',5);
    end
    if plotQTarget
        targQ0 = interp1(time,targetQ,currTime);
        linkQ0 = linkPosFun(targQ0.');
        plot(linkQ0(:,1) + xOff1_ + xOff2_,linkQ0(:,2),'color',1 - 0.4*(1-lc_),'linewidth',3,'linestyle','-'); %1 - 0.3*(1-[0,0.45,0.74])
    end
    plot(link_(:,1) + xOff1_ + xOff2_,link_(:,2),'color',lc_,'linestyle','-','linewidth',1.5);
    hipX(i+1) = link_(3,1) + xOff1_ + xOff2_;
end

plot(ha.XLim,[0,0],'k','linewidth',2)
set(gca,'xtick',hipX,'xticklabels',tVals)
xlabel('Time (s)')
% legend(lstr,'location','northwest')

end