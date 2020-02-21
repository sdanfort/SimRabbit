function [t,x,m,exitFlag] = simulateHybrid(f, S, R, x0, m0, tF, hX, hXF,...
    dim, maxstep)

if length(tF) == 2
    % tRange not final time
    t0 = tF(1);
    tF = tF(2);
else
    t0 = 0;
end

nModes = length(f);

x = cell(nModes,1);

t(1) = t0;
x{m0}(1,:) = x0';
for i = [1:m0-1, m0+1:nModes]
    x{i} = nan(1,dim(i));
end
m(1) = m0;

maxEvents = 10000;

for i = 1:maxEvents
    % Simulate trajectory in its current mode
    m_ = m(end);
    t0_ = t(end);
    x0_ = x{m_}(end,:)';
    
    % If an event is triggered at the beginning
%     IE_0_ = find(eventFun(x0_(:),S(m_,:),hX{m_},hXF{m_})>=0);    
    IE_0_ = find(eventFun(x0_(:),S(m_,:),hXF{m_})>=0);
    if ~isempty(IE_0_)
        if min(IE_0_) <= nModes 
            % Hit a guard
            mPLUS = min(IE_0_);
            xPLUS = R{m_,mPLUS}(x0_(:));

            t = [ t ; t0_];
            m = [m; mPLUS];
            x{mPLUS} = [x{mPLUS}; xPLUS'];

            for j = [1:mPLUS-1, mPLUS+1:nModes]
                x{j} = [x{j} ; nan(1,dim(j))];
            end
            continue
%         else
%             % Left the set
%             exitFlag = -1;
%             return
        end
    end
    
    if nargin < 10
        opts = odeset('events',@(~,x_) eventFun(x_(:),S(m_,:),hXF{m_}),...
            'abstol',1e-6,'reltol',1e-6);
%         opts = odeset('events',@(~,x_) eventFun(x_(:),S(m_,:),hX{m_},hXF{m_}));%,...
%                       'abstol',1e-3,'reltol',1e-2);
                  %'maxstep',1e-3);
    else
        opts = odeset('events',@(~,x_) eventFun(x_(:),S(m_,:),hXF{m_}),...
            'maxstep',maxstep);
%         opts = odeset('events',@(~,x_) eventFun(x_(:),S(m_,:),hX{m_},hXF{m_}),'maxstep',maxstep);
    end
    
    [t_,x_,TE,~,IE] = ode45(@(~,x_) f{m_}(x_(:)), [t0_,tF], x0_, opts);
    
    if ~isempty(TE)
        IE = IE(abs(TE - TE(end))<1e-4); 
    end
    
    % Append this to our trajectory
    t = [t ; t_(2:end)];
    m = [m ; m_*ones(length(t_)-1,1)];
    x{m_} = [x{m_}; x_(2:end,:)];
    for j = [1:m_-1, m_+1:nModes]
        x{j} = [x{j} ; nan(length(t_)-1,dim(j))];
    end
    
    % Determine if we left the set or hit a guard
    if (abs(t_(end) - tF)<1e-4) || isempty(IE)
        % Timed out
        exitFlag = 1;
        return  
    elseif min(IE) <= nModes 
        % Hit a guard
        mPLUS = min(IE);
        xPLUS = R{m_,mPLUS}(x_(end,:)');

        t = [ t ; t_(end)];
        m = [m; mPLUS];
        x{mPLUS} = [x{mPLUS}; xPLUS'];
        
        for j = [1:mPLUS-1, mPLUS+1:nModes]
            x{j} = [x{j} ; nan(1,dim(j))];
        end
    else
        % Left the set
        exitFlag = -1;
        return
    end
end

% Hit the maximum number of events
exitFlag = 0;
end

function [value, isterminal, direction] = eventFun(x,S,hXF)%(x,S,hX,hXF)

nModes = length(S);
eps_ = 1e-8;

S_Val = zeros(nModes,1);
for i = 1:nModes
    if isempty(S{i}(x))
        S_Val(i) = -1;
    else
        S_i = S{i}(x) + eps_;
        
        % Get minimum guard condition
        S_Val(i) = min(S_i);
%         % Use a hinge function to ensure that ALL guard conditions are met
%         negDist = min(S_i,0);
%         S_Val(i) = sum(negDist);
    end
end

% hX_Val = -(hX(x) + eps_);

hXF_Val = -ones(length(hXF),1);
if ~isempty(hXF)
    if ~iscell(hXF)
        hXF = {hXF};
    end
    for i = 1:length(hXF)
        hXF_i = hXF{i}(x) + eps_;

        % Use a hinge function to ensure that ALL guard conditions are met
        negDist = min(hXF_i,0);
        hXF_Val(i) = sum(negDist);
    end
end

value = [S_Val;
%          hX_Val;
         hXF_Val];

     
% if max(S_Val) > 100*eps_
%     disp('Event')
% end
% if max(hXF_Val) > -eps_
%     disp('Fail')
% end
% disp(max(S_Val))
isterminal = ones(length(value),1);

direction = ones(length(value),1);

end