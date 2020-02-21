classdef NdHybridModel < HybridModel
    %NdHybridModel Hybrid dynamic model with guard and reset uncertainties
    %   Detailed explanation goes here
    
    properties
        nModes  % Number of modes of the system (integer >= 1)
        nStates % Cell array with number of states for each mode
        fFN     % [nModes X 1] Cell arrays with continuous dynamics for each mode
        guFN    %                       *
        gdFN    %                       *
        hXFN    % [nModes X 1] Cell array with domain conditions for each mode
        hXFFN   % [nModes X 1] Cell array with failure conditions for each mode
        
        % Guard and reset condition for transition i->j is:
        %   S(xM,xP_NonDtr) >= 0
        %   xP_Dtr = R(xM)
        %   xP = [xP_NonDtr; xP_Dtr];
        indDtr  % [nModes X nModes] cell array of logical(nStates,1), deterministic == true
        hSFN    % [nModes X nModes] Cell array of guard conditions S(x{i},x{j}(~indExp)) 
        hRFN    % [nModes X nModes] Cell array of explicit reset maps
        
        % For simulation
        hSBarFN % hsBar(xM) > 0 if E xP s.t. hSFN(xM,xP) > 0  (e.g. max_xP(hSFN))
    end
    
    methods
        function model = HybridDynamicModel(nStates,fFN,guFN,gdFN,hXFN,hXFFN,...
                                            indDtr,hSFN,hRFN,hSBarFN)
            model.nModes = length(nStates);
            model.nStates = nStates;
            model.fFN = fFN;
            model.guFN = guFN;
            model.gdFN = gdFN;
            model.hXFN = hXFN;
            model.hXFFN = hXFFN;
            
            model.indDtr = indDtr;
            model.hSFN = hSFN;
            model.hRFN = hRFN;
            
            if ~all(indDtr)
                if (nargin < 10 || isempty(hSBarFN))
                    % Have to compute hSBar
                    % ------- To do -------
                else
                    % Don't have to compute hSBar
                    model.hSBarFN = hSBarFN;
                end
            end
        end
        
        function [t,x,m,exitFlag] = Simulate(obj,uFN,dFN,tRange,x0,m0)
            odeFN = cell(obj.nModes,1);
            for i = 1:obj.nModes
                odeFN{i} = @(x_) obj.fFN{i}(x_) + obj.guFN{i}(x_)*uFN{i}(x_) ...
                                                + obj.gDFN{i}(x_)*dFN{i}(x_);
            end
            
            if all(obj.indDtr)
                % Deterministic system
                [t,x,m,exitFlag] = simulateHybrid(odeFN, obj.SFN, obj.RFN,...
                                                  x0, m0, tRange, obj.hXFN, obj.hXFFN, obj.nStates);
            else
                % Nondeterministic system
                % ------- To do: Complete function below ------------
                [t,x,m,exitFlag] = simulateHybridNondtr(odeFN, obj.hSFN, obj.RFN, obj.hSBarFN, ...
                               x0, m0, tRange, obj.hXFN, obj.hXFFN, obj.nStates);
            end
            
        end
    end
    
end

