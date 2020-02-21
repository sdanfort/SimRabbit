function uHat = UalphaToUhatFull(uAlpha,x,xa,Tmax)
    [~, u0_, uua_] = ComputeUstarFull(x,xa,0);
    
    ua_max = min((sign(uua_)*Tmax - u0_)./uua_);
    ua_min = max((-sign(uua_)*Tmax - u0_)./uua_);
    
    if ua_max < ua_min
%         warning('No Safe Input')
        uHat = 0;
    elseif uAlpha > ua_max
        uHat = 1;
    elseif uAlpha < ua_min
        uHat = -1;
    else
        uHat = (uAlpha - (ua_max + ua_min)/2) / ((ua_max - ua_min)/2);
    end
end