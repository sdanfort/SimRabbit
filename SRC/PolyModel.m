classdef PolyModel < DynamicModel
    %POLYMODEL Dynamic model expressed with msspoly objects scaled to lie
    %between [-1,1]
    %   Allows for invariant set computation
    
    properties
        scaleParam
        dl
        
        % msspoly functions 
        x
        f
        gu
        gd
        hX
        hXF
        
        % Invariant set computation output
        invSetOut 
    end
    
    methods
        function obj = PolyModel(x,f,gu,gd,hX,hXF,scaleParam)
            % Cast everything as msspoly
            x = msspoly(x); f = msspoly(f); gu = msspoly(gu); gd = msspoly(gd); 
            hX = msspoly(hX); hXF = msspoly(hXF);
            
            % Compute matlab functions for symbolic expressions
            fFN = fn(f,x);
            guFN = fn(gu,x);
            gdFN = fn(gd,x);
            hXFN = fn(hX,x);
            hXFFN = fn(hXF,x);
            
            % Initialize dynamics model
            obj@DynamicModel(length(x),fFN,guFN,gdFN,hXFN,hXFFN);
            
            % Fill in poly objects
            obj.x = x; 
            obj.f = f; 
            obj.gu = gu;
            obj.gd = gd; 
            obj.hX = hX; 
            obj.hXF = hXF;
            
            if nargin > 6
                obj.scaleParam = scaleParam;
            end
            
            % Default moments are box moments from [-1,1]
            obj.dl = boxMoments( x, -ones(size(x)), ones(size(x)) );
        end
        
        function xUS = UnscaleState(obj,xSC)
            if isempty(obj.scaleParam)
                xUS = xSC;
            else
                switch find(size(xSC)==length(obj.x),1)
                    case 1
                        xUS = obj.scaleParam.x_US_Fun(xSC);
                    case 2
                        xUS = obj.scaleParam.x_US_Fun(xSC');
                    otherwise
                        error('State dimension mismatch')
                end
            end
        end
        
        function uUS = UnscaleController(obj,uSC)
            if isempty(obj.scaleParam)
                uUS = uSC;
            else
                switch lower(class(uSC))
                    case 'msspoly'
                        uUS = obj.scaleParam.u_US_Fun(subs(uSC,obj.x,obj.scaleParam.x_SC_Fun(obj.x)));
                    otherwise %uSC is a function
                        uUS = @(x_) obj.scaleParam.u_US_Fun(uSC(obj.scaleParam.x_SC_Fun(x_)));
                end
            end
        end
        
        function dUS = UnscaleDisturbance(obj,dSC)
            if isempty(obj.scaleParam)
                dUS = dSC;
            else
                switch lower(class(dSC))
                    case 'msspoly'
                        dUS = obj.scaleParam.d_US_Fun(subs(dSC,obj.x,obj.scaleParam.x_SC_Fun(obj.x)));
                    otherwise %uSC is a function
                        dUS = @(x_) obj.scaleParam.d_US_Fun(dSC(obj.scaleParam.x_SC_Fun(x_)));
                end
            end
        end
        
        function [obj,invSetOut] = ComputeInvariantSet(obj,options)
            % Default options
            lambda_0 = 0*obj.x;
            c_0 = 0.1;
            u_0 = zeros(size(obj.gu,2),1);
            degree = 4;
            x_TP = [];
            TP_weight = 0;
            noiseScale = 1;

            % Overwrite defaults from options object
            if nargin>1 && ~isempty(options)
                if isfield(options,'lambda_0')
                    lambda_0 = options.lambda_0;
                end
                if isfield(options,'c_0')
                    c_0 = options.c_0;
                end
                if isfield(options,'u_0')
                    u_0 = options.u_0;
                end
                if isfield(options,'degree')
                    degree = options.degree;
                end
                if isfield(options,'x_TP')
                    x_TP = options.x_TP;
                end
                if isfield(options,'TP_weight')
                    TP_weight = options.TP_weight;
                end
                if isfield(options,'noiseScale')
                    noiseScale = options.noiseScale;
                end
            end
            
            invSetOut = invariantSet(obj.x,obj.f,obj.gu,noiseScale*obj.gd,obj.hX,obj.hXF,...
                                     lambda_0,c_0,u_0,obj.dl,degree,x_TP,TP_weight);
            obj.invSetOut = invSetOut;
        end
        
        
        function [tSim,xSim] = Simulate(obj,u,d,tRange,x0)
            switch lower(class(u))
                case 'msspoly'
                    uFN = fn(u,obj.x);
                otherwise
                    uFN = u;
            end
            switch lower(class(d))
                case 'msspoly'
                    dFN = fn(d,obj.x);
                otherwise
                    dFN = d;
            end
            
            [tSim,xSim] = Simulate@DynamicModel(obj,uFN,dFN,tRange,x0);
        end
    end
    
end

