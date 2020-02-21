function RabbitModel = GenRabbitAug(dqMax)
    addpath('../../../../src');
    addpath('../autogen');
    nStates = 13;
    
    qRange = [-pi/2,pi/2; % Pitch
              -pi/2,pi/2; % Stance hip
              -pi/2,0.05;    % Stance knee
              -pi/2,pi/2  % Swing hip
              -pi/2,0.05];
    
    tauSel = [0,0,0,0;
              1,0,0,0;
              0,1,0,0;
              0,0,1,0;
              0,0,0,1];
    
    fFN = @(x) [x(6:10); % q
                MassMatrixFN(x(1:5))\CoriGravFN(x(1:5),x(6:10)); % dq
                x(12); % alpha
                0; % dalpha
                1]; % time
%     guFN = @(x) [[zeros(5,4);
%                   MassMatrixFN(x(1:5))\tauSel;
%                   zeros(3,4)],                 [zeros(11,1); 1; 0]];
    guFN = @(x) [zeros(5,5);
                  MassMatrixFN(x(1:5))\tauSel;
                  zeros(3,5)];
    hXFN = @(x) [(x(1:5) - qRange(:,1)).*(qRange(:,2) - x(1:5));
                 (x(6:10) + dqMax).*(dqMax - x(6:10))];
    hXFFN = cell(10,1);
    for i = 1:5
        hXFFN{i} = @(x)-(x(i) - qRange(i,1)).*(qRange(i,2) - x(i));
    end
    for i = 6:10
        hXFFN{i} = @(x)-(x(i) + dqMax).*(dqMax - x(i));
    end
    hXFFN{11} = @(x)[[-1,0;0,-1]*PosSwingFN(x(1:5));
                 -JYSwingFN(x(1:5))*x(6:10)];
%     hXFFN{11} = @(x) x(6) + x(7) + x(8)/2 - 0.5; % Negative theta velocity
    
    hSFN = @(x) [[1,0;0,-1]*PosSwingFN(x(1:5));
                 -JYSwingFN(x(1:5))*x(6:10)];
    RFN  = @(x) [x([1;4;5;2;3]);
                 DeltaDqFN(x(1:5))*x(6:10);
                 x(11);
                 x(12);
                 x(13)];

    RabbitModel = HybridDynamicModel({nStates},{fFN},{guFN},{[]},{hXFN},{hXFFN},{hSFN},{RFN});

end