function uAlpha = UhatToUalphaFull(uHat,x,xa,Tmax)
    [~, u0_, uua_] = ComputeUstarFull(x,xa,0);
    ua_max = min((sign(uua_)*Tmax - u0_)./uua_);
    ua_min = max((-sign(uua_)*Tmax - u0_)./uua_);
    
    if ua_max < ua_min
%         warning('No Safe Input')
        uAlpha = 0;%(ua_max + ua_min)/2;
    elseif uHat > 1
        uAlpha = ua_max;
    elseif uHat < -1
        uAlpha = ua_min;
    else
        uAlpha = uHat*(ua_max - ua_min)/2 + (ua_max + ua_min)/2;
    end
end