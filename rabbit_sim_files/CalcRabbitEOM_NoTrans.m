function CalcRabbitEOM_NoTrans
%CalcRabbitEOM Generates the equations of motion for a compass gait walker
% This MATLAB function performs the symbolic computation of the equations
% of motion for a compass gait walker with massive legs. The mass and
% length parameters of the robot are stored in the function 'CGparams.m'.
% The output of this function is the file 'CGEOM.mat', which contains
% symbolic expressions for the dynamics, guard and reset map for this
% model, as well as kinematic properties of the model for animation
% purposes.

%   More information can be found in the paper titled "Guaranteed Safe 
%   Semi-Autonomous Control of Hybrid Systems" submitted to Robotics and 
%   Automation Letters and IROS 2018.

%   Created by Nils Smit-Anseeuw (1) on 2-15-18
%   MATLAB 2017a

%   (1) Robotics and Optimization for Analysis of Human Motion
%       University of Michigan Ann Arbor
%       nilssmit@umich.edu

%% Get rabbit's parameters
param = RabbitParam;

% Gravity
g = param.g;

% Masses
MT = param.MT;
Mf = param.Mf;
Mt = param.Mt;

% Lengths
lT = param.lT;
lf = param.lf;
lt = param.lt;

% Inertias
IT = param.IT;
If = param.If;
It = param.It;

% COM Positions (From joint origin)
pMT = param.pMT;
pMf = param.pMf;
pMt = param.pMt;

% Transmission Properties
nMot = param.nMot;
IRot = param.IRot*0;

%% Define variables
syms phi alphaR betaR alphaL betaL dphi dalphaR dbetaR dalphaL dbetaL

%% Get continuous dynamics using minimal coordinates
% Minimal coordinates (right leg is on the ground)
q = [phi,alphaR,betaR,alphaL,betaL].'; 
dq = [dphi,dalphaR,dbetaR,dalphaL,dbetaL].';

%Link angles
MB_Ang = phi;
Stance_Thigh_Ang = MB_Ang + alphaR;
Stance_Shank_Ang = Stance_Thigh_Ang + betaR;
Swing_Thigh_Ang = MB_Ang + alphaL;
Swing_Shank_Ang = Swing_Thigh_Ang + betaL;

%Link positions
Stance_Foot_Pos = [0;0];
Stance_Shank_CoG = (lt - pMt)*[-sin(Stance_Shank_Ang);
                                    cos(Stance_Shank_Ang)];
Stance_Knee_Pos = lt*[-sin(Stance_Shank_Ang);
                           cos(Stance_Shank_Ang)];
Stance_Thigh_CoG = Stance_Knee_Pos +...
                       (lf - pMf)*[-sin(Stance_Thigh_Ang);
                                    cos(Stance_Thigh_Ang)];
Hip_Pos = Stance_Knee_Pos +...
              lf*[-sin(Stance_Thigh_Ang);
                   cos(Stance_Thigh_Ang)];
MB_CoG = Hip_Pos +...
             pMT*[-sin(MB_Ang);
                   cos(MB_Ang)];
Head_Pos = Hip_Pos +...
           2*pMT*[-sin(MB_Ang);
                   cos(MB_Ang)];
Swing_Thigh_CoG = Hip_Pos +...
                      pMf*[ sin(Swing_Thigh_Ang);
                           -cos(Swing_Thigh_Ang)];
Swing_Knee_Pos = Hip_Pos +...
                     lf*[ sin(Swing_Thigh_Ang);
                         -cos(Swing_Thigh_Ang)];
Swing_Shank_CoG = Swing_Knee_Pos +...
                      pMt*[ sin(Swing_Shank_Ang);
                           -cos(Swing_Shank_Ang)];
Swing_Foot_Pos = Swing_Knee_Pos +...
                     lt*[ sin(Swing_Shank_Ang);
                         -cos(Swing_Shank_Ang)];

% Link angular velocities
d_MB_Ang = jacobian(MB_Ang,q)*dq;
d_Stance_Thigh_Ang = jacobian(Stance_Thigh_Ang,q)*dq;
d_Stance_Shank_Ang = jacobian(Stance_Shank_Ang,q)*dq;
d_Swing_Thigh_Ang = jacobian(Swing_Thigh_Ang,q)*dq;
d_Swing_Shank_Ang = jacobian(Swing_Shank_Ang,q)*dq;
d_Stance_Hip_Rot = dalphaR*nMot;
d_Stance_Knee_Rot = dbetaR*nMot;
d_Swing_Hip_Rot = dalphaL*nMot;
d_Swing_Knee_Rot = dbetaL*nMot;

% Link translational velocities
d_MB_CoG = jacobian(MB_CoG,q)*dq;
d_Stance_Thigh_CoG = jacobian(Stance_Thigh_CoG,q)*dq;
d_Stance_Shank_CoG = jacobian(Stance_Shank_CoG,q)*dq;
d_Swing_Thigh_CoG = jacobian(Swing_Thigh_CoG,q)*dq;
d_Swing_Shank_CoG = jacobian(Swing_Shank_CoG,q)*dq;
%%                     
% Potential Energy
V = g*( MT * MB_CoG(2) +...
            Mf * (Stance_Thigh_CoG(2) + Swing_Thigh_CoG(2)) +...
            Mt * (Stance_Shank_CoG(2) + Swing_Shank_CoG(2)) );
V = simplify(V);

% Kinetic Energy:         
T = 0.5*( MT * sum(d_MB_CoG.^2) +...
              Mf * sum(d_Stance_Thigh_CoG.^2 + d_Swing_Thigh_CoG.^2) +...
              Mt * sum(d_Stance_Shank_CoG.^2 + d_Swing_Shank_CoG.^2) +...
              IT * d_MB_Ang^2 +...
              If * (d_Stance_Thigh_Ang^2 + d_Swing_Thigh_Ang^2) +...
              It * (d_Stance_Shank_Ang^2 + d_Swing_Shank_Ang^2) +...
              IRot * (d_Stance_Hip_Rot^2 + d_Stance_Knee_Rot^2 + ...
                      d_Swing_Hip_Rot^2 + d_Swing_Knee_Rot^2));
                  
T = simplify(T);

% Lagrangian
L = T-V;

% Partial derivatives:
dLdq   = jacobian(L,q);
dLdqdt = jacobian(L,dq);

% Compute Mass Matrix:
M = jacobian(dLdqdt.',dq);
M = simplify(M);

% Compute the coriolis and gravitational forces:
dL_dqdt_dt = jacobian(dLdqdt.',q)*dq;
f_cg = dLdq.' - dL_dqdt_dt;
f_cg = simplify(f_cg);

% fprintf('Inverting mass matrix...\n')
% 
% M_Inv = simplify(inv(M));
% 
% fprintf('Computing continuous dynamics...\n')
% 
% d_dqdt_dt_Min = (M_Inv*f_cg);

% tauSel = eye(length(q));
% tauSel = tauSel(:,2:end); % Torque on each actuated joint
% gU = simplify(M_Inv*tauSel);
% gD = gU;

% Swing foot Jacobian
J = simplify(jacobian(Swing_Foot_Pos,q.'));

%% Get discrete dynamics using angular momentum
syms dphiPlus dalphaRPlus dbetaRPlus dalphaLPlus dbetaLPlus
dqPlus = [dphiPlus,dalphaRPlus,dbetaRPlus,dalphaLPlus,dbetaLPlus].'; % Unknown, solving for these
qPlus = [phi,alphaL,betaL,alphaR,betaR].'; % Known, just swap left and right legs

hMinus = sym(zeros(5,1)); % Angular momenta before impact
hPlus =  sym(zeros(5,1));  % Angular momenta after impact

% 2d cross product
c2d = @(x1,x2) x1(1)*x2(2) - x1(2)*x2(1);

% Angular momentum conserved about impact foot
hMinus(1) = Mt*c2d(Swing_Shank_CoG  - Swing_Foot_Pos,d_Swing_Shank_CoG)  + It*d_Swing_Shank_Ang  +...
            Mf*c2d(Swing_Thigh_CoG  - Swing_Foot_Pos,d_Swing_Thigh_CoG)  + If*d_Swing_Thigh_Ang  +...
            MT*c2d(MB_CoG           - Swing_Foot_Pos,d_MB_CoG)           + IT*d_MB_Ang           +...
            Mf*c2d(Stance_Thigh_CoG - Swing_Foot_Pos,d_Stance_Thigh_CoG) + If*d_Stance_Thigh_Ang +...
            Mt*c2d(Stance_Shank_CoG - Swing_Foot_Pos,d_Stance_Shank_CoG) + It*d_Stance_Shank_Ang;
hPlus(1) = subs(...
            Mt*c2d(Swing_Shank_CoG  - Stance_Foot_Pos,d_Swing_Shank_CoG)  + It*d_Swing_Shank_Ang  +...
            Mf*c2d(Swing_Thigh_CoG  - Stance_Foot_Pos,d_Swing_Thigh_CoG)  + If*d_Swing_Thigh_Ang  +...
            MT*c2d(MB_CoG           - Stance_Foot_Pos,d_MB_CoG)           + IT*d_MB_Ang           +...
            Mf*c2d(Stance_Thigh_CoG - Stance_Foot_Pos,d_Stance_Thigh_CoG) + If*d_Stance_Thigh_Ang +...
            Mt*c2d(Stance_Shank_CoG - Stance_Foot_Pos,d_Stance_Shank_CoG) + It*d_Stance_Shank_Ang,...
            [q;dq],[qPlus;dqPlus]);
         
% Angular momentum conserved about impact knee
hMinus(2) = Mf*c2d(Swing_Thigh_CoG  - Swing_Knee_Pos,d_Swing_Thigh_CoG)  + If*d_Swing_Thigh_Ang  +...
            MT*c2d(MB_CoG           - Swing_Knee_Pos,d_MB_CoG)           + IT*d_MB_Ang           +...
            Mf*c2d(Stance_Thigh_CoG - Swing_Knee_Pos,d_Stance_Thigh_CoG) + If*d_Stance_Thigh_Ang +...
            Mt*c2d(Stance_Shank_CoG - Swing_Knee_Pos,d_Stance_Shank_CoG) + It*d_Stance_Shank_Ang;
hPlus(2) = subs(...
            Mt*c2d(Swing_Shank_CoG  - Stance_Knee_Pos,d_Swing_Shank_CoG)  + It*d_Swing_Shank_Ang  +...
            Mf*c2d(Swing_Thigh_CoG  - Stance_Knee_Pos,d_Swing_Thigh_CoG)  + If*d_Swing_Thigh_Ang  +...
            MT*c2d(MB_CoG           - Stance_Knee_Pos,d_MB_CoG)           + IT*d_MB_Ang           +...
            Mf*c2d(Stance_Thigh_CoG - Stance_Knee_Pos,d_Stance_Thigh_CoG) + If*d_Stance_Thigh_Ang,...
            [q;dq],[qPlus;dqPlus]);

% Main body angular momentum conserved
hMinus(3) = MT*c2d(MB_CoG - Hip_Pos,d_MB_CoG) + IT*d_MB_Ang;
hPlus(3)  = subs(...
            MT*c2d(MB_CoG - Hip_Pos,d_MB_CoG) + IT*d_MB_Ang,...
            [q;dq],[qPlus;dqPlus]);
        
% Angular momentum conserved about hip
hMinus(4) = Mf*c2d(Stance_Thigh_CoG - Hip_Pos,d_Stance_Thigh_CoG) + If*d_Stance_Thigh_Ang +...
            Mt*c2d(Stance_Shank_CoG - Hip_Pos,d_Stance_Shank_CoG) + It*d_Stance_Shank_Ang;
hPlus(4) = subs(...
            Mt*c2d(Swing_Shank_CoG  - Hip_Pos,d_Swing_Shank_CoG)  + It*d_Swing_Shank_Ang  +...
            Mf*c2d(Swing_Thigh_CoG  - Hip_Pos,d_Swing_Thigh_CoG)  + If*d_Swing_Thigh_Ang,...
            [q;dq],[qPlus;dqPlus]);
        
% Angular momentum conserved about liftoff knee
hMinus(5) = Mt*c2d(Stance_Shank_CoG - Stance_Knee_Pos,d_Stance_Shank_CoG) + It*d_Stance_Shank_Ang;
hPlus(5) = subs(...
            Mt*c2d(Swing_Shank_CoG  - Swing_Knee_Pos,d_Swing_Shank_CoG)  + It*d_Swing_Shank_Ang,...
            [q;dq],[qPlus;dqPlus]);

hMinus = simplify(hMinus);
hPlus = simplify(hPlus);

hMat = simplify(jacobian(hPlus,dqPlus));

fprintf('Computing impact reset...\n')
dqReset = simplify(hMat\hMinus);
DeltaDq = simplify(jacobian(dqReset,dq));

%% Save dynamics and reset map

% Link positions (for animation purposes)
linkPos = [0, 0;
           Stance_Knee_Pos.';
           Hip_Pos.';
           Head_Pos.';
           Hip_Pos.';
           Swing_Knee_Pos.';
           Swing_Foot_Pos.'];


matlabFunction(M,'vars',{q},'file','autogen/MassMatrixFN')
matlabFunction(f_cg,'vars',{q,dq},'file','autogen/CoriGravFN')
matlabFunction(Swing_Foot_Pos,'vars',{q},'file','autogen/PosSwingFN')
matlabFunction(Swing_Foot_Pos(2),'vars',{q},'file','autogen/YSwingFN')
matlabFunction(J,'vars',{q},'file','autogen/JSwingFN')
matlabFunction(J(2,:),'vars',{q},'file','autogen/JYSwingFN')
matlabFunction(DeltaDq,'vars',{q},'file','autogen/DeltaDqFN')
matlabFunction(linkPos,'vars',{q},'file','autogen/LinkPosFN')

%% Save analytic representations in a file

save('RabbitEOM_Full','q','dq','M','f_cg','Swing_Foot_Pos','J','DeltaDq','linkPos')
end
