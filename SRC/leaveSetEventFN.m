function [val,isterm,dir] = leaveSetEventFN(hXFN,xVal)
    val = hXFN(xVal);
    isterm = ones(length(val),1);
    dir = -ones(length(val),1);
end