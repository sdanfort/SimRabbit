function uMaskFN = maskInputInterpHybrid(uSafe, vMask, x, vThresh)
    nModes = length(vMask);
    
    uMaskFN = cell(nModes,1);

    vSpline = spline([vThresh/2,vThresh],[0 0 1 0]);

    for i = 1:nModes
        if ~isempty(uSafe{i})
            vMaskFN = fn(vMask{i},x{i});
            uSafeFN = fn(uSafe{i},x{i});

            uMaskFN{i} = @(x_,u_) switchFN(vMaskFN(x_), vThresh, vSpline, u_, uSafeFN(x_));
        end
    end
end

function [u_,v,z_] = switchFN(v,thresh,spline,u_0,u_1)
    if v < thresh/2
        z_ = 0;
    elseif v < thresh
        z_ = ppval(spline,v);
    else
        z_ = 1;
    end 
    
    u_ = max(-1,min(1,z_*u_0 + (1-z_)*u_1));
end