function [RdTheta,ERdTheta] = FitPolyReset(x,deg,nSamp,nTest,xRange)
eps_ = 0.0001;


%% Sample dynamics
thMinusUS_ = xRange(1,2);
[dthUS_,aUS_,daUS_] = ndgrid(linspace(xRange(2,1),xRange(2,2),nSamp),...
                       linspace(xRange(3,1),xRange(3,2),nSamp),...
                       linspace(xRange(4,1),xRange(4,2),nSamp));
xUS_ = [thMinusUS_+0*dthUS_(:),dthUS_(:),aUS_(:),daUS_(:)];

[dThPlusInitUS_] = sampleReset(xUS_);

%% Scale state to fit in 1-by-1 box
xMid_ = mean(xRange,2).';
xScale_ = diff(xRange,[],2).'/2;

xInit_ = (xUS_ - repmat(xMid_,length(xUS_),1))./repmat(xScale_,length(xUS_),1);

dThPlusInit_ = (dThPlusInitUS_ - xMid_(2))/xScale_(2);

%% Generate poly fits
% Define variables

monoms = monomials(x(2:4),0:deg);

monomsInit_ = full(msubs(monoms,x,xInit_.'));

[dthTest_,aTest_,daTest_] = ndgrid(linspace(-1,1,nTest),linspace(-1,1,nTest),linspace(-1,1,nTest));
xTest_ = [1+0*dthTest_(:),dthTest_(:),aTest_(:),daTest_(:)];
nTest = size(xTest_,1);
monomsTest_ = full(msubs(monoms,x,xTest_.'));

xTestUS_ = xTest_.*repmat(xScale_,nTest,1) + repmat(xMid_,nTest,1);
[dThPlusTestUS_] = sampleReset(xTestUS_);

dThPlusTest_ = (dThPlusTestUS_ - xMid_(2))/xScale_(2);

bm_ = boxMoments(x,-ones(4,1),ones(4,1));
int_ = double(bm_(monoms));

%%      Get poly fit for dThetaPlus
fprintf('Fitting dThetaPlus\n')

dThPlusLB = AdaptiveSampledPolyFit(monoms,int_,monomsInit_,dThPlusInit_,monomsTest_,dThPlusTest_,-1,eps_);
dThPlusUB = AdaptiveSampledPolyFit(monoms,int_,monomsInit_,dThPlusInit_,monomsTest_,dThPlusTest_,1,eps_);

%%      Plot the results

[dthPlot_,aPlot_] = meshgrid(linspace(-1,1,50));
daPlot_ = linspace(-1,1,3);
nPlot = length(dthPlot_(:));

for i = 1:length(daPlot_)
    xPlot_ = [1+0*dthPlot_(:),dthPlot_(:),aPlot_(:),daPlot_(i)*ones(nPlot,1)];
    xPlotUS_ = xPlot_.*repmat(xScale_,nPlot,1) + repmat(xMid_,nPlot,1);
    
    dThPlusPlotUS_ = sampleReset(xPlotUS_);
    
    dThPlusPlot_ = reshape((dThPlusPlotUS_ - xMid_(2))/xScale_(2),size(dthPlot_));
        
    LBPlot_ = reshape(double(msubs(dThPlusLB,x,xPlot_.')).',size(dthPlot_));
    UBPlot_ = reshape(double(msubs(dThPlusUB,x,xPlot_.')).',size(dthPlot_));  
    
    figure(i+5),clf,hold on
    surf(dthPlot_,aPlot_,dThPlusPlot_,'facecolor','b')
    surf(dthPlot_,aPlot_,LBPlot_,'facecolor','r')
    surf(dthPlot_,aPlot_,UBPlot_,'facecolor','g')
        
    view(3)
end

%% Export polynomials
RdTheta = (dThPlusUB + dThPlusLB)/2;
ERdTheta = (dThPlusUB - dThPlusLB)/2;

end
