function [ thSim, xSim, dxSim, x0 ] = load_gait_library( idx )

%load the gait library and return an x-trajectory.

%Input is idx, between 1 and 46.
%Lower idx-values correspond to narrower step widtths.

%load FROST gaits
obj_ = load('Data/RabbitGaits_og');
rabbitGaits = obj_.rabbitGaits(5:end,:);

nTrajs = size( rabbitGaits, 1 );
nStates = 5;
nNodes = length(rabbitGaits{1,3}(1).tspan);

qTraj = zeros(nTrajs,nNodes,nStates);
dqTraj = zeros(nTrajs,nNodes,nStates);
uTraj = zeros(nTrajs, nNodes,4);

timeMeshgrid = zeros(nTrajs,nNodes);

%note: we don't actually have to load all of them if you just want one 
%for this function output hahaha
%but here's how you could if you needed
for i = 1:nTrajs

    %getting into the 5-DOF form we have depicted in the paper
    qTarget_ = rabbitGaits{i,3}(1).states.x( 3:end, : ).';
    qTarget_( :,1 ) = -qTarget_(:,1);
    qTarget_( :,[2,4] ) = pi - qTarget_(:,[2,4]);
    qTarget_( :,[3,5] ) = -qTarget_(:,[3,5]);
    
        %adjust the foot height a bit so it's actually hitting the ground
    footHeight_ = @(phiMod_) YSwingFN(qTarget_(1,:).' +...
        [phiMod_;0;0;0;0]) - YSwingFN(qTarget_(end,:).' + [phiMod_; 0;...
        0; 0; 0]);
    
    phiMod_ = fsolve( footHeight_, 0, optimoptions( 'fsolve', 'Display',...
        'none' ) );

    qTarget_ = qTarget_ + [ phiMod_, 0, 0, 0, 0 ];
    
    %define theta, the stance leg angle of robot
    th_tmp = -qTarget_(:,1) - qTarget_(:,2) - qTarget_(:,3)/2;
    
    qTraj(i,:,:) = qTarget_;
    
    %note: there is currently no dq traj or input. I'm adding here:
    dqTarget_ = rabbitGaits{i,3}(1).states.dx( 3:end, : ).';
    dqTarget_ = -dqTarget_;
    dqTraj(i,:,:) = dqTarget_;
    
    uTarget_ = rabbitGaits{i,3}(1).inputs.u.';
    uTraj(i,:,:) = -uTarget_;
    %we also probably have to reverse direction of these
    
    thetaTraj(i,:) = th_tmp;
    timeMeshgrid(i,:) = round(rabbitGaits{i,3}(1).tspan.'*1000)/1000;
    
end

%define tSim and xSim straight from FROST data
xSim = [ qTraj( idx, :, 1 )', qTraj( idx, :, 2 )', qTraj( idx, :, 3 )'...
    qTraj( idx, :, 4 )', qTraj( idx, :, 5 )' ];
dxSim = [ dqTraj( idx, :, 1 )', dqTraj( idx, :, 2 )', dqTraj( idx, :, 3 )'...
    dqTraj( idx, :, 4 )', dqTraj( idx, :, 5 )' ];
thSim = thetaTraj( idx, : )';
x0 = [ qTraj( idx, 1, 1 )', qTraj( idx, 1, 2 )', qTraj( idx, 1, 3 )'...
    qTraj( idx, 1, 4 )', qTraj( idx, 1, 5 )', dqTraj( idx, 1, 1 )',...
    dqTraj( idx, 1, 2 )', dqTraj( idx, 1, 3 )'...
    dqTraj( idx, 1, 4 )', dqTraj( idx, 1, 5 )' ];