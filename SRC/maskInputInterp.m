function uMaskFN = maskInputInterp(uSafe, vMask, x, vThresh)
    vMaskFN = fn(vMask,x);
    uSafeFN = fn(uSafe,x);
    
    vSpline = spline([vThresh/2,vThresh],[0 0 1 0]);
    
%     zFN = @(x_) 2/(1+exp(-vMaskFN(x_)/vThresh)) - 1;

    uMaskFN = @(x_,u_) switchFN(vMaskFN(x_), vThresh, vSpline, u_, uSafeFN(x_));
end

function u_ = switchFN(val,thresh,spline,u_0,u_1)
    if val < thresh/2
        z_ = 0;
    elseif val < thresh
        z_ = ppval(spline,val);
    else
        z_ = 1;
    end 
    
    u_ = z_*u_0 + (1-z_)*u_1;
end