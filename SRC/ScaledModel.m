classdef ScaledModel
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        UnscaledDynamicModel
        ScaledDynamicModel
        ScaleParams
    end
    
    methods
        function out = ScaledModel(USMod, SCMod, SCParams)
            out.UnscaledDynamicModel = USMod;
            out.ScaledDynamicModel = SCMod;
            out.ScaleParams = SCParams;
        end
        
        function [t,xUS] = simulateScaled(Model,x0SC,tRange)
            
        end
        
        function [val,xUS] = evaluate(Model,fun,xSC)
            
        end
    end
    
end

