classdef DynamicModel
    %DYNAMICMODEL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        nStates
        fFN
        guFN
        gdFN
        hXFN
        hXFFN
    end
    
    methods
        function model = DynamicModel(nStates,fFN,guFN,gdFN,hXFN,hXFFN)
            model.nStates = nStates;
            model.fFN = fFN;
            model.guFN = guFN;
            model.gdFN = gdFN;
            if nargin > 4
                model.hXFN = hXFN;
                model.hXFFN = hXFFN;
            end
        end
        
        function [tSim,xSim] = Simulate(obj,uFN,dFN,tRange,x0)
            opts = odeset('events',@(~,x_) leaveSetEventFN(obj.hXFN,x_));

            [tSim,xSim] = ode45(@(t_,x_) obj.fFN(x_) +...
                                         obj.guFN(x_)*max(-1,min(1,uFN(x_))) +...
                                         obj.gdFN(x_)*max(-1,min(1,dFN(x_))),...
                                         tRange,x0,opts);
        end
    end
    
end

