function [v_0, v_0_LowDeg, uLQR, uLQR_LowDeg] = RabbitLQR(polyMod,xPer,mPer,uPer,plotFlag)
if ~exist('plotFlag','var')
    plotFlag = 0;
end

    nModes = max(mPer);
thInd = 1;
nonThInd = [2,3,4];

%% Fit polynomials to the periodic trajectory
xPerUS = mergeModes(polyMod.UnscaleState(xPer,mPer),mPer);
uPerUS = zeros(length(xPerUS),1);
for i = 1:nModes
    uPerUS(mPer==i) = polyMod.scaleParam.u_US_Fun{i}(uPer{i}(mPer==i));
end
figure,plot(xPerUS(:,thInd),xPerUS(:,nonThInd)),hold on
plot(xPerUS(:,thInd),uPerUS)

fitOrder = 3;
x_0 = cell(nModes,1);
x_0_Lin = x_0;
u_0 = cell(nModes,1);
u_0_Lin = u_0;
theta = sym('th');
monom = theta.^[fitOrder:-1:0].';

for i = 1:nModes
    x_0{i} = sym(zeros(3,1));
    
    thetaPer = xPer{i}(mPer == i,thInd);
    xNonThetaPer = xPer{i}(mPer == i,nonThInd);
    xNonThetaFit_ = zeros(length(thetaPer),3);
    xNonThetaLin_ = zeros(length(thetaPer),3);
    
    for j = 1:3
        lsqC_ = thetaPer.^(fitOrder:-1:0);
        lsqd_ = xNonThetaPer(:,j);
        lsqB_ = [lsqC_(1,:);lsqC_(end,:)];
        lsqc_ = [lsqd_(1);lsqd_(end)];
        coeff_ = lsqlin(lsqC_,lsqd_,[],[],lsqB_,lsqc_).';
        x_0{i}(j) = coeff_*monom;
        
        xNonThetaFit_(:,j) = polyval(coeff_,thetaPer);
        
%         coeffLin_ = polyfit(thetaPer,xNonThetaPer(:,j),1);
        x_0_Lin{i}(j) = xNonThetaPer(1,j) + (theta - thetaPer(1))/(thetaPer(end)-thetaPer(1))*(xNonThetaPer(end,j)-xNonThetaPer(1,j));
        xNonThetaLin_(:,j) = subs(x_0_Lin{i}(j),theta,thetaPer);
    end
    
    xFit_ = zeros(length(thetaPer),length(nonThInd)+1);
    xFit_(:,thInd) = thetaPer;
    xFit_(:,nonThInd) = xNonThetaFit_;
    xFitUS_ = polyMod.UnscaleState(xFit_,i);
    
    xLin_ = zeros(length(thetaPer),length(nonThInd)+1);
    xLin_(:,thInd) = thetaPer;
    xLin_(:,nonThInd) = xNonThetaLin_;
    xLinUS_ = polyMod.UnscaleState(xLin_,i);
    
    plot(xFitUS_(:,thInd),xFitUS_(:,nonThInd),'--k')
    plot(xLinUS_(:,thInd),xLinUS_(:,nonThInd),':k')
    
    uP_ = uPer{i}(mPer == i);
    lsqC_ = thetaPer.^(fitOrder:-1:0);
    lsqd_ = uP_;
    lsqB_ = [lsqC_(1,:);lsqC_(end,:)];
    lsqc_ = [lsqd_(1);lsqd_(end)];
    coeff_ = lsqlin(lsqC_,lsqd_,[],[],lsqB_,lsqc_).';
%     coeff_ = polyfit(thetaPer,uPer{i}(mPer == i).',fitOrder);
    u_0{i} = coeff_*monom;
    
    u_0_Lin{i} = uP_(1) + (theta - thetaPer(1))/(thetaPer(end)-thetaPer(1))*uP_(end);
    
    uFit_ = polyval(coeff_,thetaPer);
    uFitUS_ = polyMod.scaleParam.u_US_Fun{i}(uFit_);
    
    plot(xFitUS_(:,thInd),uFitUS_,'--k')
end

%% Get time-varying linearized dynamics
A = cell(nModes,1);
Afn = cell(nModes,1);
B = cell(nModes,1);
Bfn = cell(nModes,1);

xHat = sym('xh',[3,1]);

for i = 1:nModes
    xSub = sym(zeros(4,1));
    xSub(nonThInd) = xHat + x_0{i};
    xSub(thInd) = theta;
    
    fSub_ = polyMod.fFN{i}(xSub);
    guSub_ = polyMod.guFN{i}(xSub);
    dth_ = fSub_(thInd);
    
    fHat = (fSub_(nonThInd) + guSub_(nonThInd,:)*u_0{i})/dth_ - jacobian(x_0{i},theta);
    
    gHat = guSub_(nonThInd,:)/dth_;
    
    A{i} = subs(jacobian(fHat,xHat),xHat,zeros(3,1));
    Afn{i} = matlabFunction(A{i},'vars',theta);
    B{i} = subs(gHat,xHat,zeros(3,1));
    Bfn{i} = matlabFunction(B{i},'vars',theta);
end

%% Simulate t-v dynamics and compare to full order
xLTV = repcell({[xHat;theta]},nModes,1);
fLTV = cell(nModes,1);
guLTV = cell(nModes,1);
gdLTV = repcell({[]},nModes,1);
hXLTV = repcell({[]},nModes,1);
hXFLTV = repcell({[]},nModes,1);
hSLTV = cell(nModes,nModes);
RLTV = cell(nModes,nModes);

thetaLim = linspace(-1,1,nModes+1);
for i = 1:nModes
    fLTV{i} = [A{i}*xHat;
                       1];
    guLTV{i} = [B{i};
                   0];
               
    if i < nModes
        hSLTV{i,i+1} = theta - thetaLim(i+1);
        RLTV{i,i+1} = [xHat;
                       theta];
    end
end

ltvMod = SymHybridModel(xLTV,fLTV,guLTV,gdLTV,hXLTV,hXFLTV,hSLTV,RLTV);

uLTV = repcell({-3*xHat(2) - 1*xHat(3)},nModes,1);
x0LTV = [0.0;0.1;0;-1];
[tOutLTV,xOutLTV,mOutLTV] = ltvMod.Simulate(uLTV,repcell({[]},nModes,1),[0,2],x0LTV,1);
xOutLTV = mergeModes(xOutLTV,mOutLTV);
figure,hold on,plot(tOutLTV-1,xOutLTV);

uPoly =cell(nModes,1);
xHatPoly = cell(nModes,1);

for i = 1:nModes
    u_0f_ = matlabFunction(u_0{i},'vars',theta);
    x_0f_ = matlabFunction(x_0{i},'vars',theta);
    uLTVf_ = matlabFunction(uLTV{i},'vars',{xHat});
    thetaPoly_ = polyMod.x{i}(thInd);
    xHatPoly{i} = polyMod.x{i}(nonThInd) - x_0f_(thetaPoly_);
    uPoly{i} = u_0f_(thetaPoly_) + uLTVf_(xHatPoly{i});
end
x0Poly = zeros(4,1);
x0Poly(nonThInd) = x0LTV(1:3) + subs(x_0{1},theta,-1);
x0Poly(thInd) = -1;

[~,xOutPoly,mOutPoly] = polyMod.Simulate(uPoly,repcell({@(t)0},nModes,1),[0,0.4],x0Poly,1);
for i = 1:max(mOutPoly)
    plot(xOutPoly{i}(mOutPoly == i,thInd),xOutPoly{i}(mOutPoly == i,nonThInd) - double(subs(x_0{i},theta,xOutPoly{i}(mOutPoly == i,thInd).')).','--k')
end

%% Compute LQR regulator by integrating Riccati backwards
% Q = 0*eye(3);    % State cost
% R = 1000;         % Control cost
% F = 100*eye(3); % Terminal cost
% 
Q = eye(3);    % State cost
R = 5;         % Control cost
F = 100*eye(3); % Terminal cost

P = cell(nModes,1);
K = cell(nModes,1);
nDisc = 15;
thetaSteps = cell(nModes,1);

for i = nModes:-1:1
    thetaSteps{i} = linspace(thetaLim(i),thetaLim(i+1),nDisc);
    h = diff(thetaSteps{i}(1:2));
    P{i} = zeros(3,3,nDisc);
    K{i} = zeros(1,3,nDisc);
    if i == nModes
        P{i}(:,:,end) = F;
    else
        P{i}(:,:,end) = P{i+1}(:,:,1);
    end
    
    for j = nDisc:-1:1
        A_ = Afn{i}(thetaSteps{i}(j));
        B_ = Bfn{i}(thetaSteps{i}(j));
        P_ = P{i}(:,:,j);
        
        K{i}(:,:,j) = R^-1*(B_.'*P_);
        
        if j > 1
            dP_ = A_.'*P_ + P_*A_ - P_*B_*R^-1*B_.'*P_ + Q;

            P{i}(:,:,j-1) = P_ + h*dP_;
        end
    end
end

%% Visualize the 1 level set of P(theta)
if plotFlag
    figure, hold on
    camlight; lighting phong
    xSlice_ = 0;
    for i = 1:nModes
        thSt_ = thetaSteps{i}(1:end);
        [thPlot,xPlot1,xPlot2] = ndgrid(thSt_,linspace(-1,1,101),linspace(-1,1,101));
        V_ = 0*xPlot1;
        for j = 1:length(thSt_)
            inds = (thPlot == thSt_(j));
            x_ = [xPlot1(inds),xPlot2(inds),xSlice_ + 0*xPlot2(inds)];
            P_ = P{i}(:,:,j);
            V_(inds) = sum((x_*P_).*x_,2);
        end

        patch_ = patch(isosurface(thPlot,xPlot1,xPlot2,V_,1));

        patch_.FaceColor = [240,159,30]/256;
        patch_.EdgeColor = 'none';

%         isonormals(thPlot,xPlot1,xPlot2,V_,patch_);
    end
%     
%     dphi_ = 0;
%     for i = 1:nModes
%         thSt_ = thetaSteps(1:end);
%         [phi_,th_,dth_] = ndgrid(thSt_);
%         V_ = 0*phi_;
%         for j = 1:length(th_)
%             inds = (th_ == thSt_(j));
%             x_ = [phi_(inds),zeros(sum(inds(:)),1)+dphi_,dth_(inds)];
%             P_ = P{i}(:,:,(j-1)+1);
%             V_(inds) = sum((x_*P_).*x_,2);
%         end
% 
%         patch_ = patch(isosurface((th_-1+2*i)/(2*nModes),phi_,dth_,V_,1));
% 
%         patch_.FaceColor = [240,159,30]/256;
%         patch_.EdgeColor = 'none';
% 
%         isonormals((th_-1+2*i)/(2*nModes),phi_,dth_,V_,patch_);
%     end
    
    axis equal
end
%% Interpolate P to get v_0
highDeg = 4;
lowDeg = 2;
monom = polyMod.x{1}(thInd).^[highDeg:-1:0].';
monomLowDeg = polyMod.x{1}(thInd).^[lowDeg:-1:0].';

v_0 = cell(nModes,1);
v_0_LowDeg = v_0;

uLQR = cell(nModes,1);
uLQR_LowDeg = uLQR;

PfitFun = cell(nModes,1);

% lsqMat = (1-thetaSteps(2:end-1).^2).'*(1-thetaSteps(2:end-1).^2);

for i = 1:nModes
    Pfit_ = msspoly(zeros(3,3));
    PfitLowDeg_ = msspoly(zeros(3,3));
    Kfit_ = msspoly(zeros(1,3));
    KfitLowDeg_ = msspoly(zeros(1,3));
    for j = 1:3
        for k = 1:3
            lsqC_ = thetaSteps{i}.'.^(highDeg:-1:0);
            lsqd_ = squeeze(P{i}(j,k,:));
            lsqB_ = [lsqC_(1,:);lsqC_(end,:)];
            lsqc_ = [lsqd_(1);lsqd_(end)];
            coeff_ = lsqlin(lsqC_,lsqd_,[],[],lsqB_,lsqc_).';
            Pfit_(j,k) = coeff_*monom;
            
            % For low degree, ensure fit matches at boundaries
            lsqC_ = thetaSteps{i}.'.^(lowDeg:-1:0);
            lsqd_ = squeeze(P{i}(j,k,:));
            lsqB_ = [lsqC_(1,:);lsqC_(end,:)];
            lsqc_ = [lsqd_(1);lsqd_(end)];
            coeff_ = lsqlin(lsqC_,lsqd_,[],[],lsqB_,lsqc_).';
            PfitLowDeg_(j,k) = coeff_*monomLowDeg;
        end
        
        lsqC_ = thetaSteps{i}.'.^(highDeg:-1:0);
        lsqd_ = squeeze(K{i}(1,j,:));
        lsqB_ = [lsqC_(1,:);lsqC_(end,:)];
        lsqc_ = [lsqd_(1);lsqd_(end)];
        coeff_ = lsqlin(lsqC_,lsqd_,[],[],lsqB_,lsqc_).';
%         coeff_ = polyfit(thetaSteps.',squeeze(K{i}(1,j,:)),highDeg);
        Kfit_(1,j) = coeff_*monom;
        
        lsqC_ = thetaSteps{i}.'.^(lowDeg:-1:0);
        lsqd_ = squeeze(K{i}(1,j,:));
        lsqB_ = [lsqC_(1,:);lsqC_(end,:)];
        lsqc_ = [lsqd_(1);lsqd_(end)];
        coeff_ = lsqlin(lsqC_,lsqd_,[],[],lsqB_,lsqc_).';
%         coeff_ = polyfit(thetaSteps.',squeeze(K{i}(1,j,:)),lowDeg);
        KfitLowDeg_(1,j) = coeff_*monomLowDeg;
    end
    
    PfitFun{i} = fn(Pfit_,polyMod.x{1}(thInd));
    
    v_0{i} = 1 - xHatPoly{i}.'*Pfit_*xHatPoly{i};
    
    u_0fn_ = matlabFunction(u_0{i},'vars',theta);
    
    uLQR{i} = u_0fn_(polyMod.x{i}(thInd)) - Kfit_*xHatPoly{i};
    
    % Get low order v
    x_0Linf_ = matlabFunction(x_0_Lin{i},'vars',theta);
    xHatLinPoly = polyMod.x{i}(nonThInd) - x_0Linf_(polyMod.x{i}(thInd)).';
    
    v_0_LowDeg{i} = 1 - xHatLinPoly.'*PfitLowDeg_*xHatLinPoly;
    
    u_0_Linf_ = matlabFunction(u_0_Lin{i},'vars',theta);
    uLQR_LowDeg{i} = u_0_Linf_(polyMod.x{i}(thInd)) - KfitLowDeg_*xHatLinPoly;
end

%% Get max u values
for i = 1:nModes
    uMax_ = polyMax(uLQR_LowDeg{i},polyMod.x{i},polyMod.hX{i},8);
    uMin_ = -polyMax(-uLQR_LowDeg{i},polyMod.x{i},polyMod.hX{i},8);
    fprintf('Mode %i: uMax = %4.3f, uMin = %4.3f\n',i,uMax_,uMin_)
end

%% Plot comparison to P
if plotFlag
    figure, hold on
    camlight; lighting phong
    xSlice_ = 0;
    for i = 1:nModes
        thSt_ = thetaSteps{i}(1:end);
        [thPlot,xPlot1,xPlot2] = ndgrid(thSt_,linspace(-1,1,101),linspace(-1,1,101));
        V_ = 0*xPlot1;
        for j = 1:length(thSt_)
            inds = (thPlot == thSt_(j));
            x_ = [xPlot1(inds),xPlot2(inds),xSlice_ + 0*xPlot2(inds)];
            P_ = PfitFun{i}(thSt_(j));
            V_(inds) = sum((x_*P_).*x_,2);
        end

        patch_ = patch(isosurface(thPlot,xPlot1,xPlot2,V_,1));

        patch_.FaceColor = [240,159,30]/256;
        patch_.EdgeColor = 'none';

%         isonormals(thPlot,xPlot1,xPlot2,V_,patch_);
    end
    
    axis equal
    

%% Plot in full space
    figure, hold on
    camlight; lighting phong
    
    xSlice_ = 0;
    xSym = sym('x',[4,1]);
    
    for i = 1:nModes
        thSt_ = thetaSteps{i}(1:end);
        [thPlot,dthPlot,aPlot] = ndgrid(thSt_,linspace(-1,1,101),linspace(-1,1,101));
        vfn_ = fn(v_0{i},polyMod.x{i});
        vfn_ = matlabFunction(vfn_(xSym),'vars',{xSym});
        v_ = vfn_([thPlot(:),dthPlot(:),aPlot(:),0*aPlot(:)+xSlice_].');
        v_ = reshape(v_,size(thPlot));

        patch_ = patch(isosurface(thPlot,dthPlot,aPlot,v_,0));

        patch_.FaceColor = [240,159,30]/256;
        patch_.EdgeColor = 'none';

%         isonormals(thPlot,dthPlot,aPlot,v_,patch_);
    end
    
    axis equal
    

%% Plot low degree
    figure, hold on
    camlight; lighting phong
    
    xSlice_ = 0;
    xSym = sym('x',[4,1]);
    
    for i = 1:nModes
        thSt_ = thetaSteps{i}(1:end);
        [thPlot,dthPlot,aPlot] = ndgrid(thSt_,linspace(-1,1,101),linspace(-1,1,101));
        vfn_ = fn(v_0_LowDeg{i},polyMod.x{i});
        vfn_ = matlabFunction(vfn_(xSym),'vars',{xSym});
        v_ = vfn_([thPlot(:),dthPlot(:),aPlot(:),0*aPlot(:)+xSlice_].');
        v_ = reshape(v_,size(thPlot));

        patch_ = patch(isosurface(thPlot,dthPlot,aPlot,v_,0));

        patch_.FaceColor = [240,159,30]/256;
        patch_.EdgeColor = 'none';

%         isonormals(thPlot,dthPlot,aPlot,v_,patch_);
    end
    
    axis equal
    
end

end
