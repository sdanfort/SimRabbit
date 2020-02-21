function CalcRabbitEOM
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

syms xStanceFoot yStanceFoot dxStanceFoot dyStanceFoot

%% Get continuous dynamics using minimal coordinates
% Minimal coordinates (right leg is on the ground)
q = [xStanceFoot,yStanceFoot,phi,alphaR,betaR,alphaL,betaL].'; 
dq = [dxStanceFoot,dyStanceFoot,dphi,dalphaR,dbetaR,dalphaL,dbetaL].';

%Link angles
MB_Ang = phi;
Stance_Thigh_Ang = MB_Ang + alphaR;
Stance_Shank_Ang = Stance_Thigh_Ang + betaR;
Swing_Thigh_Ang = MB_Ang + alphaL;
Swing_Shank_Ang = Swing_Thigh_Ang + betaL;

%Link positions
Stance_Foot_Pos = [xStanceFoot;yStanceFoot];
Stance_Shank_CoG = Stance_Foot_Pos +...
                       (lt - pMt)*[-sin(Stance_Shank_Ang);
                                    cos(Stance_Shank_Ang)];
Stance_Knee_Pos = Stance_Foot_Pos +...
                       lt*[-sin(Stance_Shank_Ang);
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

% Link angular velocities (assuming stance foot pos is fixed)
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

%% Get continuous dynamics using minimal coordinates (swing foot position = 0)
q_Min = [phi,alphaR,betaR,alphaL,betaL].';
dq_Min = [dphi,dalphaR,dbetaR,dalphaL,dbetaL].';

L_Min = subs(L,[xStanceFoot,yStanceFoot,dxStanceFoot,dyStanceFoot],[0,0,0,0]);

% Partial derivatives:
dLdq_Min   = jacobian(L_Min,q_Min);
dLdqdt_Min = jacobian(L_Min,dq_Min);

% Compute Mass Matrix:
M_Min = jacobian(dLdqdt_Min.',dq_Min);
M_Min = simplify(M_Min);

% Compute the coriolis and gravitational forces:
dL_dqdt_dt_Min = jacobian(dLdqdt_Min.',q_Min)*dq_Min;
f_cg_Min = dLdq_Min.' - dL_dqdt_dt_Min;
f_cg_Min = simplify(f_cg_Min);

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
Swing_Foot_Pos_Min = subs(Swing_Foot_Pos,[xStanceFoot,yStanceFoot,dxStanceFoot,dyStanceFoot],[0,0,0,0]);
J_Min = simplify(jacobian(Swing_Foot_Pos_Min,q_Min.'));

%% Get discrete dynamics using floating body dynamics
fprintf('Calculating Reset...\n')
% Partial derivatives:
dLdqdt = jacobian(L,dq);

% Compute Mass Matrix:
M = jacobian(dLdqdt.',dq);
fprintf('   Simplifying M\n')
M = simplify(M);

fprintf('   Simplifying J\n')
J = simplify(jacobian(Swing_Foot_Pos,q.'));

fprintf('   Computing M_tilde\n')
M_Tilde = (J*(M\(J.')));

% Update map in velocity (floating body)
fprintf('   Computing DeltaDq\n')
DeltaDq = (eye(length(q)) - M\(J.'*(M_Tilde\J)));

DeltaDq_Min = [1,0,0,0,0;
               0,0,0,1,0;
               0,0,0,0,1;
               0,1,0,0,0;
               0,0,1,0,0]*DeltaDq(3:end,3:end);

DeltaDq_Min = subs(DeltaDq_Min,[xStanceFoot,yStanceFoot],[0,0]);
%% Save dynamics and reset map
sprintf('Exporting...\n')
% Link positions (for animation purposes)
linkPos_Min = subs([0, 0;
               Stance_Knee_Pos.';
               Hip_Pos.';
               Head_Pos.';
               Hip_Pos.';
               Swing_Knee_Pos.';
               Swing_Foot_Pos.'],[xStanceFoot,yStanceFoot,dxStanceFoot,dyStanceFoot],[0,0,0,0]);


matlabFunction(M_Min,'vars',{q_Min},'file','autogen/MassMatrixFN')
matlabFunction(f_cg_Min,'vars',{q_Min,dq_Min},'file','autogen/CoriGravFN')
matlabFunction(Swing_Foot_Pos_Min,'vars',{q_Min},'file','autogen/PosSwingFN')
matlabFunction(Swing_Foot_Pos_Min(2),'vars',{q_Min},'file','autogen/YSwingFN')
matlabFunction(J_Min,'vars',{q_Min},'file','autogen/JSwingFN')
matlabFunction(J_Min(2,:),'vars',{q_Min},'file','autogen/JYSwingFN')
matlabFunction(DeltaDq_Min,'vars',{q_Min},'file','autogen/DeltaDqFN')
matlabFunction(linkPos_Min,'vars',{q_Min},'file','autogen/LinkPosFN')

%% Save analytic representations in a file

%save('RabbitEOM_Full','q','dq','M','f_cg','Swing_Foot_Pos','J','DeltaDq','linkPos')
end
