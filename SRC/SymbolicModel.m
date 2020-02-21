classdef SymbolicModel < DynamicModel
    %SYMBOLICMODEL Dynamic model expressed with symbolic objects
    %   Allows for scaling and tayloring of dynamics
    
    properties
        % Symbolic expressions for each dynamic element
        x
        f
        gu
        gd
        hX
        hXF
    end
    
    methods
        function obj = SymbolicModel(x,f,gu,gd,hX,hXF)
            if nargin < 5
                hX = [];
            end
            if nargin < 6
                hXF = [];
            end
            
            % Cast everything as sym
            x = sym(x); f = sym(f); gu = sym(gu); gd = sym(gd); 
            hX = sym(hX); hXF = sym(hXF);
            
            % Compute matlab functions for symbolic expressions
            fFN = matlabFunction(f, 'vars', {x});
            guFN = matlabFunction(gu, 'vars', {x});
            gdFN = matlabFunction(gd, 'vars', {x});
            hXFN = matlabFunction(hX, 'vars', {x});
            hXFFN = matlabFunction(hXF, 'vars', {x});
            
            % Initialize dynamics model
            obj@DynamicModel(length(x),fFN,guFN,gdFN,hXFN,hXFFN);
            
            % Fill in symbolic objects
            obj.x = x; 
            obj.f = f; 
            obj.gu = gu;
            obj.gd = gd; 
            obj.hX = hX; 
            obj.hXF = hXF;
        end
        
        function modelSC = ScaleTaylor(obj, xRange, uRange, dRange, degree)
            x_US_Sym = obj.x;
            fUS = obj.f;
            guUS = obj.gu;
            gdUS = obj.gd;

            x_SC = msspoly('x',length(x_US_Sym));
            [fSC, guSC, gdSC, scaleParam] = scaleTaylorDynamics(x_US_Sym, x_SC, fUS, guUS, gdUS, xRange, uRange, dRange, degree);

            if ~isempty(obj.hX)
                hXSC = subs(obj.hX,x_US_Sym,scaleParam.x_US_Expr);
            else
                hXSC = (x_SC + ones(size(x_SC))).*(ones(size(x_SC)) - x_SC);
            end
            if ~isempty(obj.hXF)
                hXFSC = subs(obj.hXF,x_US_Sym,scaleParam.x_US_Expr);
            else
                hXFSC = [];
            end

            modelSC = PolyModel(x_SC,fSC,guSC,gdSC,hXSC,hXFSC,scaleParam);
        end
        
        function [tSim,xSim] = Simulate(obj,u,d,tRange,x0)
            switch lower(class(u))
                case 'sym'
                    uFN = matlabFunction(u, 'vars', {obj.x});
                otherwise
                    uFN = u;
            end
            switch lower(class(d))
                case 'sym'
                    dFN = matlabFunction(d, 'vars', {obj.x});
                otherwise
                    dFN = d;
            end
            
            [tSim,xSim] = Simulate@DynamicModel(obj,uFN,dFN,tRange,x0);
        end
    end
    
end

