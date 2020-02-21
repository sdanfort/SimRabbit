function [val,isterm,dir] = leaveSetEvent(hX,x,xVal)
    val = double(subs(hX,x,xVal));
    isterm = ones(length(hX),1);
    dir = -ones(length(hX),1);
end