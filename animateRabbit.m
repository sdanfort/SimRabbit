function animateRabbit(time,states,linkPosFun )

%pause value for drawing the links
n = .1;

%line color for plotting
lc_ = [0,0,0];

time = time-time(1);

%detect gait events (heel strike!)
eventInds = find(~diff(time))+1;
eventInds = [1;eventInds(:);size(states,1)];

%get our first link position
link0 = linkPosFun(states(1,:).');

figure(1);
clf;
hold on;
xlim([-0.5 0.5]);
hl = plot( link0(:,1), link0(:,2), 'color', lc_, 'linewidth', 2);
ha = gca;
axis equal

hold on;
plot( ha.XLim, [0,0], 'k' )

yl = ha.YLim;
xl = ha.XLim;
drawnow;
% pause(n);

tic
xOffset = link0(1,1);

for i = 1:length(eventInds)-1
    currTime = toc;
    inds = eventInds(i):(eventInds(i+1)-1);
    while currTime < time(inds(end))
        
       st_ = interp1(time(inds),states(inds,:),currTime);
       link_ = linkPosFun(st_.');
       
       set(hl,'xdata',link_(:,1)+xOffset);
       set(hl,'ydata',link_(:,2));
       set(ha,'ylim',yl);
       set(ha,'xlim',link_(3,1)+xOffset+xl)
       drawnow
       pause(n);
       currTime = toc;
       
    end
    
    %update the offset so we translate our x-direction
    link_ = linkPosFun(states(inds(end),:).');
    xOffset = xOffset + link_(end,1);
    
end

end
    