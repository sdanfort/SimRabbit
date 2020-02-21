function dl = RabbitPPMoments(x,dThMax,nModes)
if nargin < 3
    nModes = 1;
end

if ~iscell(dThMax)
    dThMoment = repcell({@(p) subs(integral(p,x(2)),x(2),dThMax) -...
        subs(integral(p,x(2)),x(2),-1)},nModes,1);
else
    dThMoment = cell(nModes,1);
    for i = 1:length(dThMax)
        dThMoment{i} = @(p) subs(integral(p,x(2)),x(2),dThMax{i}) -...
            subs(integral(p,x(2)),x(2),-1);
    end
end
if nModes == 1
    m_ = boxMoments( x([1,3,4]), -ones(3,1), ones(3,1) );
    dl = @(p) m_(dThMoment{1}(p));
else
    thetaLim = linspace(-1,1,nModes+1);
    dl = cell(nModes,1);
    for i = 1:nModes
        m_ = boxMoments( x([1,3,4]), [thetaLim(i);-ones(2,1)],...
            [thetaLim(i+1);ones(2,1)] );
        dl{i} = @(p) m_(dThMoment{i}(p));
    end
end    
    
end