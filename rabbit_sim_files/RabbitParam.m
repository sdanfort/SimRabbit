function param = RabbitParam
%% Get the compass gait parameters
% Gravity
param.g = 9.81;

% Masses
param.MT = 12;
param.Mf = 6.8;
param.Mt = 3.2;

% Lengths
param.lT = 0.63;
param.lf = 0.4;
param.lt = 0.4;

% Inertias
param.IT = 1.33;
param.If = 0.47;
param.It = 0.2;

% COM Positions (From joint origin)
param.pMT = 0.24;
param.pMf = 0.11;
param.pMt = 0.24;

% Transmission Properties
param.nMot = 50; % Motor gear ratio
param.IRot = 0.000332; % Rotor inertia