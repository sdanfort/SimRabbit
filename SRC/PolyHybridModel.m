classdef PolyHybridModel < HybridDynamicModel
    %POLYHYBRIDMODEL Hybrid dynamic model expressed with msspoly objects scaled to lie
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
        hS
        R
        
        % Invariant set computation output
        invSetOut 
    end
    
    methods
        function obj = PolyHybridModel(x,f,gu,gd,hX,hXF,hS,R,scaleParam)
            nModes = length(x);
            
            nStates =  zeros(nModes,1);
            fFN = cell(nModes,1); guFN = cell(nModes,1); gdFN =...
                cell(nModes,1);
            hXFN = cell(nModes,1); hXFFN = cell(nModes,1); dl =...
                cell(nModes,1);
            hSFN = cell(nModes,nModes); RFN = cell(nModes,nModes);
            
            for i = 1:nModes
                nStates(i) = length(x{i});
                
                % Cast everything as msspoly
                x{i} = msspoly(x{i}); f{i} = msspoly(f{i}); gu{i} =...
                    msspoly(gu{i}); 
                gd{i} = msspoly(gd{i}); hX{i} = msspoly(hX{i}); 
                for j = 1:length(hXF{i})
                    hXF{i}{j} = msspoly(hXF{i}{j});
                end

                % Compute matlab functions for msspoly expressions
                fFN{i} = fn(f{i},x{i});
                if ~isempty(gu{i})
                    guFN{i} = fn(gu{i},x{i});
                else
                    guFN{i} = [];
                end
                if ~isempty(gd{i})
                    gdFN{i} = fn(gd{i},x{i});
                else
                    gdFN{i} = [];
                end
                hXFN{i} = fn(hX{i},x{i});
                for j = 1:length(hXF{i})
                    hXFFN{i}{j} = fn(hXF{i}{j},x{i});
                end
                
                for j = 1:nModes
                    hS{i,j} = msspoly(hS{i,j}); R{i,j} = msspoly(R{i,j});
                    hSFN{i,j} = fn(hS{i,j},x{i});
                    RFN{i,j} = fn(R{i,j},x{i});
                end
            end
            
            % Initialize dynamics model
            obj@HybridDynamicModel(nStates,fFN,guFN,gdFN,hXFN,hXFFN,...
                hSFN,RFN);
            
            % Fill in poly objects
            obj.x = x; 
            obj.f = f; 
            obj.gu = gu;
            obj.gd = gd; 
            obj.hX = hX; 
            obj.hXF = hXF;
            obj.hS = hS;
            obj.R = R;
            
            if nargin > 8
                obj.scaleParam = scaleParam;
            end
            
            % Default moments are box moments from [-1,1]
            for i = 1:nModes
                obj.dl{i} = boxMoments( x{i}, -ones(size(x{i})),...
                    ones(size(x{i})));
            end
        end
    
        function xSC = ScaleState(obj,xUS,m)
            if ~iscell(xUS)
                flip = 0;
                if size(xUS,1) ~= 4 && size(xUS,2) == 4
                    xUS = xUS';
                    flip = 1;
                end
                
                xSC = obj.scaleParam.x_SC_Fun{m}(xUS);
                
                if flip
                    xSC = xSC';
                end
            else
                xSC = cell(obj.nModes,1);
                
                for i = 1:obj.nModes
                    xUS_ = xUS{i};
                    xSC_Fun_ = obj.scaleParam.x_SC_Fun{i};
                    
                    flip = 0;
                    if size(xUS_,1) ~= 4 && size(xUS_,2) == 4
                        xUS_ = xUS_';
                        flip = 1;
                    end
                    
                    xSC{i} = nan(size(xUS_));
                    xSC{i}(:,m==i) = xSC_Fun_(xUS_(:,m==i));
                    
                    if flip
                        xSC{i} = xSC{i}';
                    end
                end
            end
        end
        
        function xUS = UnscaleState(obj,xSC,m)
            if ~iscell(xSC)
                flip = 0;
                if size(xSC,1) ~= 4 && size(xSC,2) == 4
                    xSC = xSC';
                    flip = 1;
                end
                
                xUS = obj.scaleParam.x_US_Fun{m}(xSC);
                
                if flip
                    xUS = xUS';
                end
            else
                xUS = cell(obj.nModes,1);
                
                for i = 1:obj.nModes
                    xUS_ = xSC{i};
                    xUS_Fun_ = obj.scaleParam.x_US_Fun{i};
                    
                    flip = 0;
                    if size(xUS_,1) ~= 4 && size(xUS_,2) == 4
                        xUS_ = xUS_';
                        flip = 1;
                    end
                    
                    xUS{i} = nan(size(xUS_));
                    xUS{i}(:,m==i) = xUS_Fun_(xUS_(:,m==i));
                    
                    if flip
                        xUS{i} = xUS{i}';
                    end
                end
            end
        end
        
        function uSC = ScaleController(obj,uUS)
            uSC = cell(obj.nModes,1);
            
            for i = 1:obj.nModes
                uUS_ = uUS{i};
                uSC_Fun_ = obj.scaleParam.u_SC_Fun{i};
                xUS_Expr_ = obj.scaleParam.x_US_Expr_Poly{i};
                switch lower(class(uUS_))
                    case 'msspoly'
                        uSC{i} = uSC_Fun_(subs(uUS_,obj.x{i},xUS_Expr_));
                    case 'function_handle'
                        uSC{i} = uSC_Fun_(uUS_(xUS_Expr_));
                end
            end
        end
        
        function uUS = UnscaleController(obj,uSC)
            uUS = cell(obj.nModes,1);
            
            for i = 1:obj.nModes
                uSC_ = uSC{i};
                uUS_Fun_ = obj.scaleParam.u_US_Fun{i};
                xSC_Expr_ = obj.scaleParam.x_SC_Fun{i}(obj.x);
                switch lower(class(uSC_))
                    case 'msspoly'
                        uUS{i} = uUS_Fun_(subs(uSC_,obj.x,xSC_Expr_));
                    case 'function_handle'
                        uUS{i} = uUS_Fun_(uSC_(xSC_Expr_));
                end
            end
        end
        
        function [tSim,xSim,mSim,exitFlag] = Simulate(obj,u,d,tRange,x0,m0)
            uFN = cell(obj.nModes,1);
            dFN = cell(obj.nModes,1);

            for i = 1:obj.nModes
                switch lower(class(u{i}))
                    case 'msspoly'
                        uFN{i} = fn(u{i},obj.x{i});
                    otherwise
                        uFN{i} = u{i};
                end
                switch lower(class(d{i}))
                    case 'msspoly'
                        dFN{i} = fn(d{i},obj.x{i});
                    otherwise
                        dFN{i} = d{i};
                end
            end

            [tSim,xSim,mSim,exitFlag] = Simulate@HybridDynamicModel(obj,...
                uFN,dFN,tRange,x0,m0);
        end
        
        function obj = BoundModelDifference(obj,fullModel,degGD)
            nModes = obj.nModes;
            minCoeff = 1e-7;
            
            gdNew = obj.gd;

            for i = 1:nModes
                nStates = obj.nStates(i);
                
                for j = 1:nStates
                    gD_j = boundPoly(obj.f{i}(j) - fullModel.f{i}(j),...
                        obj.x{i}, obj.hX{i}, obj.dl{i}, degGD, minCoeff);
                    
                    for k = 1:size(obj.gu{i},2)
                            gD_j = gD_j + boundPoly(obj.gu{i}(j) -...
                                fullModel.gu{i}(j), obj.x{i}, obj.hX{i},...
                                obj.dl{i}, degGD, minCoeff);
                    end
                    
                    gD_max = polyMax(gD_j,obj.x{i},obj.hX{i},degGD);
                    
                    fprintf('Max difference in mode %i, state %i: %4.3f\n',...
                        i,j,gD_max)
                    
                    if gD_max >= 1e-3
                        if isempty(gdNew{i})
                            gdNew{i} = msspoly(zeros(nStates,1));
                        else
                            gdNew{i} = [gdNew{i},msspoly(zeros(nStates,1))];
                        end
                        gdNew{i}(j,end) = gD_j;
                    end
                end
            end
            
            obj.gd = gdNew;
        end
        
        function [obj,invSetOut] = ComputeInvariantSet(obj, options)
            nModes = length(obj.x);
            if nargin < 2
                options = struct;
            end
            
            % Default options
            lambda_0 = repcell(msspoly(0),nModes,1);
            sigma_0 = repcell(msspoly(1),nModes,nModes);
            c_0 = 0.1;
            u_0 = cell(nModes,1);
            for i = 1:nModes
                u_0{i} = msspoly(zeros(size(obj.gu{i},2),1));
            end
            degree = 4;
            x_TP = repcell({[]},nModes,1);
            TP_weight = 0;
            x_Con = repcell({[]},nModes,1);
            noiseScale = 1;
            plotFun = [];
            vFirst = 0;
            v_0 = repcell({[]},nModes,1);
            sExp = repcell({[]},nModes,nModes);
            bang = 0;
            convIter = 1;
            inputScale = 1;
            
            % Overwrite defaults from options object
            if nargin > 1
                optNames = fieldnames(options);
                for i = 1:length(optNames)
                    n = optNames{i};
                    switch n
                        case 'lambda_0'
                            lambda_0 = options.lambda_0;
                        case 'sigma_0'
                            sigma_0 = options.sigma_0;
                        case 'c_0'
                            c_0 = options.c_0;
                        case 'u_0'
                            u_0 = options.u_0;
                        case 'degree'
                            degree = options.degree;
                        case 'x_TP'
                            x_TP = options.x_TP;
                        case 'TP_weight'
                            TP_weight = options.TP_weight;
                        case 'x_Con'
                            x_Con = options.x_Con;
                        case 'noiseScale'
                            noiseScale = options.noiseScale;
                        case 'plotFun'
                            plotFun = options.plotFun;
                        case 'v_0'
                            v_0 = options.v_0;
                            vFirst = 1;
                        case 'sExp'
                            sExp = options.sExp;
                        case 'bang'
                            bang = options.bang;
                        case 'convIter'
                            convIter = options.convIter;
                        case 'inputScale'
                            inputScale = options.inputScale;
                        otherwise
                            warning(['No recognized option: ' n])
                    end
                end
            end
            
            % Scale the noise
            gd_ = cell(nModes,1);
            if noiseScale > 0
                for i = 1:nModes
                    gd_{i} = obj.gd{i}*noiseScale;
                end
            end
            
            if bang
                invSetOut = invariantSetHybridBang2Q(obj.x,obj.f,...
                    obj.gu,gd_,obj.hS,obj.R,obj.hX,obj.hXF,lambda_0,...
                    sigma_0,c_0,u_0,obj.dl,degree,x_TP,TP_weight,...
                    x_Con,plotFun,sExp);
            else
            
                if vFirst
                    invSetOut = invariantSetHybridV0(obj.x,obj.f,...
                        obj.gu,gd_,obj.hS,obj.R,obj.hX,obj.hXF,v_0,...
                        sigma_0,c_0,obj.dl,degree,x_TP,TP_weight,...
                        x_Con,plotFun);
                else
                    invSetOut = invariantSetHybrid(obj.x,obj.f,...
                        obj.gu,gd_,obj.hS,obj.R,obj.hX,obj.hXF,...
                        lambda_0,sigma_0,c_0,u_0,obj.dl,degree,...
                        x_TP,TP_weight,x_Con,plotFun,sExp,...
                        convIter,inputScale);
                end
            end
            obj.invSetOut = invSetOut;
        end
    end
end

