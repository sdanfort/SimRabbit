function RabbitModel = GenRabbit(dqMax,fixReset)
    if nargin < 2
        fixReset = 1;
    end

%     addpath('../../../src');
%     addpath('autogen');
    nStates = 10;
    
    qRange = [-pi/2,pi/2; % Pitch
              -pi/2,pi/2; % Swing hip
              -pi/2,0;    % Swing knee
              -pi/2,pi/2  % Stance hip
              -pi/2,0];
    
    tauSel = [0,0,0,0;
              1,0,0,0;
              0,1,0,0;
              0,0,1,0;
              0,0,0,1];
    
    fFN = @(x) [x(6:10);
                MassMatrixFN(x(1:5))\CoriGravFN(x(1:5),x(6:10))];
    guFN = @(x) [zeros(5,4);
                 MassMatrixFN(x(1:5))\tauSel];
    hXFN = @(x) [(x(1:5) - qRange(:,1)).*(qRange(:,2) - x(1:5));
                 (x(6:10) + dqMax).*(dqMax - x(6:10))];
    hSFN = @(x) [[1,0;0,-1]*PosSwingFN(x(1:5) + 0.0013);
                 -JYSwingFN(x(1:5))*x(6:10)];
    RFN  = @(x) [x([1;4;5;2;3]);
                 DeltaDqFN(x(1:5))*x(6:10)];

    RabbitModel = HybridDynamicModel({nStates},{fFN},{guFN},{[]},{hXFN},{[]},{hSFN},{RFN});

end