classdef SymHybridModel < HybridDynamicModel
    %SYMHYBRIDMODEL Hybrid dynamic model expressed with symbolic objects
    %   Allows for invariant set computation
    
    properties
        % Symbolic functions 
        x
        f
        gu
        gd
        hX
        hXF
        hS
        R
        
        % Invariant set computation output
        invSetOut 
    end
    
    methods
        function obj = SymHybridModel(x,f,gu,gd,hX,hXF,hS,R)
            nModes = length(x);
            
            nStates =  zeros(nModes,1);
            fFN = cell(nModes,1); guFN = cell(nModes,1); gdFN = cell(nModes,1);
            hXFN = cell(nModes,1); hXFFN = cell(nModes,1); dl = cell(nModes,1);
            hSFN = cell(nModes,nModes); RFN = cell(nModes,nModes);
            
            for i = 1:nModes
                nStates(i) = length(x{i});
                
                % Cast everything as msspoly
                x{i} = sym(x{i}); f{i} = sym(f{i}); gu{i} = sym(gu{i}); 
                gd{i} = sym(gd{i}); hX{i} = sym(hX{i}); 
                for j = 1:length(hXF{i})
                    hXF{i}{j} = sym(hXF{i}{j});
                end

                % Compute matlab functions for symbolic expressions
                fFN{i} = matlabFunction(f{i},'vars',x(i));
                if ~isempty(gu{i})
                    guFN{i} = matlabFunction(gu{i},'vars',x(i));
                else
                    guFN{i} = [];
                end
                if ~isempty(gd{i})
                    gdFN{i} = matlabFunction(gd{i},'vars',x(i));
                else
                    gdFN{i} = [];
                end
                hXFN{i} = matlabFunction(hX{i},'vars',x(i));
                for j = 1:length(hXF{i})
                    hXFFN{i}{j} = matlabFunction(hXF{i}{j},'vars',x(i));
                end
                
                % Default moments are box moments from [-1,1]
                dl{i} =  boxMoments( x{i}, -ones(size(x{i})), ones(size(x{i})) );
                
                for j = 1:nModes
                    hS{i,j} = sym(hS{i,j}); R{i,j} = sym(R{i,j});
                    hSFN{i,j} = matlabFunction(hS{i,j},'vars',x(i));
                    RFN{i,j} = matlabFunction(R{i,j},'vars',x(i));
                end
            end
            
            % Initialize dynamics model
            obj@HybridDynamicModel(nStates,fFN,guFN,gdFN,hXFN,hXFFN,hSFN,RFN);
            
            % Fill in symbolic objects
            obj.x = x; 
            obj.f = f; 
            obj.gu = gu;
            obj.gd = gd; 
            obj.hX = hX; 
            obj.hXF = hXF;
            obj.hS = hS;
            obj.R = R;
        end
        
        function modelSC = ScaleTaylor(obj, xRange, uRange, dRange, degree, x0Taylor)
            xUS = obj.x;
            fUS = obj.f;
            guUS = obj.gu;
            gdUS = obj.gd;
            hSUS = obj.hS;
            RUS = obj.R;

            xSC = cell(obj.nModes,1);
            for i = 1:obj.nModes
                xSC{i} = msspoly('x',length(xUS{i}));
            end
            
            if ~exist('x0Taylor','var')
                x0Taylor = cell(obj.nModes,1);
                for m = 1:obj.nModes
                    x0Taylor{m} = mean(xRange{m},2);
                end
            end
                
            [fSC, guSC, gdSC, hSSC, RSC, scaleParam] = scaleTaylorDynamicsHybrid(xUS, xSC, fUS, guUS, gdUS, hSUS, RUS, xRange, uRange, dRange, degree, x0Taylor);
            
            %% Scale domain (if necessary)
            if ~isstruct(degree) || ~isfield(degree,'hX')
                deghX = 2;
            else
                deghX = degree.hX;
            end
            
            hXSC = cell(obj.nModes,1);
            for m = 1:obj.nModes
                if ~isempty(obj.hX{m})
                    hXSC_Sym = subs(obj.hX{m},xUS{m},scaleParam.x_US_Expr{m});
                    hXSC{m} = msspoly(zeros(size(hXSC_Sym)));
                    
                    for i = 1:length(hXSC{m}(:))
                        hXSC{m}(i) = msstaylor(hXSC_Sym(i),scaleParam.x_SC_Sym{m},xSC{m},x0Taylor{m},deghX+1);
                    end
                else
                    hXSC{m} = (xSC{m} + ones(size(xSC{m}))).*(ones(size(xSC{m})) - xSC{m});
                end
            end
               
            %% Scale failure set (if necessary)
            if ~isstruct(degree) || ~isfield(degree,'hXF')
                deghXF = 2;
            else
                deghXF = degree.hXF;
            end
            
            hXFSC = cell(obj.nModes,1);
            for m = 1:obj.nModes
                if ~isempty(obj.hXF{m})
                    for j = 1:length(obj.hXF{m})
                        hXFSC_Sym = subs(obj.hXF{m}{j},xUS{m},scaleParam.x_US_Expr{m});
                        hXFSC{m}{j} = msspoly(zeros(size(hXFSC_Sym)));

                        for k = 1:length(hXFSC{m}{j}(:))
                            hXFSC{m}{j}(k) = msstaylor(hXFSC_Sym(k),scaleParam.x_SC_Sym{m},xSC{m},x0Taylor{m},deghXF+1);
                        end
                    end
                else
                    hXFSC{m} = [];
                end
            end
            
            modelSC = PolyHybridModel(xSC,fSC,guSC,gdSC,hXSC,hXFSC,hSSC,RSC,scaleParam);
        end
        
        function xUS = UnscaleState(obj,xSC)
            % ------------------------ To Do ------------------------------
%             if isempty(obj.scaleParam)
%                 xUS = xSC;
%             else
%                 switch find(size(xSC)==length(obj.x),1)
%                     case 1
%                         xUS = obj.scaleParam.x_US_Fun(xSC);
%                     case 2
%                         xUS = obj.scaleParam.x_US_Fun(xSC');
%                     otherwise
%                         error('State dimension mismatch')
%                 end
%             end
        end
        
        function uUS = UnscaleController(obj,uSC)
            % ------------------------ To Do ------------------------------
%             if isempty(obj.scaleParam)
%                 uUS = uSC;
%             else
%                 switch lower(class(uSC))
%                     case 'msspoly'
%                         uUS = obj.scaleParam.u_US_Fun(subs(uSC,obj.x,obj.scaleParam.x_SC_Fun(obj.x)));
%                     otherwise %uSC is a function
%                         uUS = @(x_) obj.scaleParam.u_US_Fun(uSC(obj.scaleParam.x_SC_Fun(x_)));
%                 end
%             end
        end
        
        function dUS = UnscaleDisturbance(obj,dSC)
            % ------------------------ To Do ------------------------------
%             if isempty(obj.scaleParam)
%                 dUS = dSC;
%             else
%                 switch lower(class(dSC))
%                     case 'msspoly'
%                         dUS = obj.scaleParam.d_US_Fun(subs(dSC,obj.x,obj.scaleParam.x_SC_Fun(obj.x)));
%                     otherwise %uSC is a function
%                         dUS = @(x_) obj.scaleParam.d_US_Fun(dSC(obj.scaleParam.x_SC_Fun(x_)));
%                 end
%             end
        end
        
        function [obj,invSetOut] = ComputeInvariantSet(obj,options)
            % ------------------------ To Do ------------------------------
            % Default options
%             lambda_0 = 0*obj.x;
%             c_0 = 0.1;
%             u_0 = 0*obj.x;
%             degree = 4;
%             x_TP = [];
%             TP_weight = 0;
%             noiseScale = 1;
% 
%             % Overwrite defaults from options object
%             if isfield(options,'lambda_0')
%                 lambda_0 = options.lambda_0;
%             end
%             if isfield(options,'c_0')
%                 c_0 = options.c_0;
%             end
%             if isfield(options,'u_0')
%                 u_0 = options.u_0;
%             end
%             if isfield(options,'degree')
%                 degree = options.degree;
%             end
%             if isfield(options,'x_TP')
%                 x_TP = options.x_TP;
%             end
%             if isfield(options,'TP_weight')
%                 TP_weight = options.TP_weight;
%             end
%             if isfield(options,'noiseScale')
%                 noiseScale = options.noiseScale;
%             end
%             
%             invSetOut = invariantSet(obj.x,obj.f,obj.gu,noiseScale*obj.gd,obj.hX,obj.hXF,...
%                                      lambda_0,c_0,u_0,obj.dl,degree,x_TP,TP_weight);
%             obj.invSetOut = invSetOut;
        end
        
        
        function [tSim,xSim,mSim,exitFlag] = Simulate(obj,u,d,tRange,x0,m0)
            uFN = cell(obj.nModes,1);
            dFN = cell(obj.nModes,1);
            
            for i = 1:obj.nModes
                switch lower(class(u{i}))
                    case 'sym'
                        uFN{i} = matlabFunction(u{i}, 'vars', {obj.x{i}});
                    otherwise
                        uFN{i} = u{i};
                end
                switch lower(class(d{i}))
                    case 'sym'
                        dFN{i} = matlabFunction(d{i}, 'vars', {obj.x{i}});
                    otherwise
                        dFN{i} = d{i};
                end
            end
            
            [tSim,xSim,mSim,exitFlag] = Simulate@HybridDynamicModel(obj,uFN,dFN,tRange,x0,m0);
        end
    end
    
end

