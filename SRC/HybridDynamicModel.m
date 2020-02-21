classdef HybridDynamicModel
    %HYBRIDDYNAMICMODEL Hybrid dynamic model with no guard or reset uncertainties
    %   Detailed explanation goes here
    
    properties
        nModes  % Number of modes of the system (integer >= 1)
        nStates % Array with number of states for each mode
        fFN     % [nModes X 1] Cell arrays with continuous dynamics for each mode
        guFN    %                       *
        gdFN    %                       *
        hXFN    % [nModes X 1] Cell array with domain conditions for each mode
        hXFFN   % [nModes X 1] Cell array with failure conditions for each mode
        
        % Guard and reset condition for transition i->j is:
        %   S{i,j}(xM) >= 0
        %   xP = R{i,j}(xM)
        hSFN    % [nModes X nModes] Cell array of guard conditions S(x{i},x{j})) 
        RFN    % [nModes X nModes] Cell array of explicit reset maps
    end
    
    methods
        function model = HybridDynamicModel(nStates,fFN,guFN,gdFN,hXFN,hXFFN,...
                                            hSFN,RFN)
            model.nModes = length(nStates);
            model.nStates = nStates;
            model.fFN = fFN;
            model.guFN = guFN;
            model.gdFN = gdFN;
            model.hXFN = hXFN;
            model.hXFFN = hXFFN;
            
            model.hSFN = hSFN;
            model.RFN = RFN;
        end
        
        function [t,x,m,exitFlag] = Simulate(obj,uFN,dFN,tRange,x0,m0)
            odeFN = cell(obj.nModes,1);
            for i = 1:obj.nModes
                if isempty(obj.guFN{i})
                    uTerm = @(x_) 0;
                else
                    uTerm = @(x_) obj.guFN{i}(x_)*uFN{i}(x_);
                end
                if isempty(obj.gdFN{i})
                    dTerm = @(x_) 0;
                else
                    dTerm = @(x_) obj.gdFN{i}(x_)*dFN{i}(x_);
                end
                odeFN{i} = @(x_) obj.fFN{i}(x_) + uTerm(x_) + dTerm(x_);
            end
            
            % Deterministic system
            [t,x,m,exitFlag] = simulateHybrid(odeFN, obj.hSFN, obj.RFN,...
                                              x0, m0, tRange, obj.hXFN,...
                                              obj.hXFFN, obj.nStates);
            
        end
    end
    
end

