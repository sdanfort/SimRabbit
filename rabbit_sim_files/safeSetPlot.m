classdef safeSetPlot
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        vFN
        hAx
        hCV
        hCF
        thVSC
        dthVSC
        thFUS
        dthFUS
        nModes
        xSize
        xMid
        Tmax
        eps_v
        plotFeatures
    end
    
    methods
        function obj = safeSetPlot(polyMod,nTh,ndTh,dThMax,eps_v,Tmax,hAx,plotFeatures)
            if nargin < 4
                dThMax = 3;
            end
            if nargin < 5
                eps_v = 0.1;
            end
            obj.eps_v = eps_v;
            if nargin < 6
                Tmax = 50;
            end
            obj.Tmax = Tmax;
            if nargin < 7
                figure(7),clf
                obj.hAx = gca; hold on
            else
                obj.hAx = hAx; hold on
            end
            if nargin < 8
                plotFeatures = [1,1,1];
            end
            colormap([0,0.5,0.01; % V boundary
                      1- 0.5*(1-[0,0.5,0.01]); % V interior
                      [0.94, 0.33, 0.31]]) % FL infeasible
            caxis([0,2*eps_v])
            
            nModes = polyMod.nModes;
            xSize = polyMod.scaleParam.xSize;
            xMid = polyMod.scaleParam.xMid;
            obj.nModes = nModes; obj.xSize = xSize{1}; obj.xMid = xMid{1};
            thetaLim = linspace(-1,1,nModes+1);
            
            %% Get sample points in theta and dtheta
            thVSC = cell(nModes,1);
            dthVSC = thVSC;
            
            nThPerMode = round(nTh/nModes);
            thVal_ = linspace(thetaLim(1),thetaLim(2),nThPerMode);
            dthVal = linspace(-1,1,ndTh);
            [thVSC{1},dthVSC{1}] = meshgrid(thVal_,dthVal);
            thVUS = xSize{1}(1)*thVSC{1} + xMid{1}(1);
            dthVUS = xSize{1}(2)*dthVSC{1} + xMid{1}(2);
            
            for i = 2:nModes
                thVal_ = linspace(thetaLim(i),thetaLim(i+1),nThPerMode);
                thVal_ = thVal_(2:end);
                [thVSC{i},dthVSC{i}] = meshgrid(thVal_,dthVal);
                thVUS = [thVUS, xSize{i}(1)*thVSC{i} + xMid{i}(1)];
                dthVUS = [dthVUS,xSize{i}(2)*dthVSC{i} + xMid{i}(2)];
            end
            obj.thVSC = thVSC;
            obj.dthVSC = dthVSC;
            
            
            [obj.thFUS,obj.dthFUS] = meshgrid(linspace(xMid{1}(1) - xSize{1}(1),xMid{nModes}(1) + xSize{nModes}(1),nTh),linspace(1,dThMax,ndTh));
            
            %% Get the v function
            if plotFeatures(1)
                xSym = sym('x',[4,1]);
                v = polyMod.invSetOut.v;
                vFN = cell(nModes,1);
                for i = 1:nModes
                    vfn_ = fn(v{i},polyMod.x{i});
                    vFN{i} = matlabFunction(vfn_(xSym),'vars',{xSym});
                end
                obj.vFN = vFN;

                %% Plot the v contour
                v_ = obj.vEval(0,0);
                [~,hCV] = contourf(thVUS,dthVUS,v_,[0,eps_v]);
                obj.hCV = hCV;
            end
            
            %% Plot the fv contour
            if plotFeatures(2)
                patch([xMid{1}(1) - xSize{1}(1),xMid{nModes}(1) + xSize{nModes}(1),xMid{nModes}(1) + xSize{nModes}(1),xMid{1}(1) - xSize{1}(1)],...
                      [-0.39,-0.39,0,0],[0.94, 0.33, 0.31]);
            end
            if plotFeatures(3)
                fv_ = obj.fvEval(0,0)+2*eps_v;
                [~,obj.hCF] = contourf(obj.thFUS,obj.dthFUS,fv_,2*eps_v*[1,1]);
            end
            
            %% Add labels
            xlabel('$$\theta$$ (rad)','interpreter','latex')
            ylabel('$$\dot{\theta}$$ (rad/s)','interpreter','latex')
            
            obj.plotFeatures = plotFeatures;
        end
        
        function update(obj,aUS,daUS)
            % Update the v contour
            if obj.plotFeatures(1)
                aSC = (aUS-obj.xMid(3))/obj.xSize(3);
                daSC = (daUS-obj.xMid(4))/obj.xSize(4);
                vVal = obj.vEval(aSC,daSC);

                obj.hCV.ZData = vVal;
            end
            % Update the f contour
            if obj.plotFeatures(3)
                fv_ = obj.fvEval(aUS,daUS)+2*obj.eps_v;
                obj.hCF.ZData = fv_;
            end
        end
        
        function v_ = vEval(obj, aSC, daSC)
            v_ = [];
            for i = 1:obj.nModes
                xSC_ = [obj.thVSC{i}(:),obj.dthVSC{i}(:),0*obj.thVSC{i}(:)+aSC,0*obj.thVSC{i}(:)+daSC];
                v_ = [v_,reshape(obj.vFN{i}(xSC_.'),size(obj.dthVSC{i}))];
            end
        end
        
        function fv_ = fvEval(obj, aUS, daUS)
            xHatVals = [obj.thFUS(:),obj.dthFUS(:),0*obj.thFUS(:)+aUS,0*obj.thFUS(:)+daUS];
            failVals = sampleDynamics(xHatVals,obj.Tmax,0);
            fv_ = reshape(failVals,size(obj.thFUS));
        end
    end
end