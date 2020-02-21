function RabbitModel = GenRabbitAugPer(dqMax)
    addpath('../../../src');
    addpath('autogen');
    nStates = 13;
    
    qRange = [-pi/3,pi/3; % Pitch
              -pi,pi; % Stance hip
              -pi/2,0.05;    % Stance knee
              -pi,pi  % Swing hip
              -pi/2,0.05];
    
    tauSel = [0,0,0,0;
              1,0,0,0;
              0,1,0,0;
              0,0,1,0;
              0,0,0,1];
    
    fFN = @(x) [x(6:10);
                MassMatrixFN(x(1:5))\CoriGravFN(x(1:5),x(6:10));
                x(12);
                0;
                1];
    guFN = @(x) [[zeros(5,4);
                  MassMatrixFN(x(1:5))\tauSel;
                  zeros(3,4)],                 [zeros(11,1); 1; 0]];
    hXFN = @(x) [(x(1:5) - qRange(:,1)).*(qRange(:,2) - x(1:5));
                 (x(6:10) + dqMax).*(dqMax - x(6:10))];
    hXFFN = cell(10,1);
    for i = 1:5
        hXFFN{i} = @(x)-(x(i) - qRange(i,1)).*(qRange(i,2) - x(i));
    end
    for i = 6:10
        hXFFN{i} = @(x)-(x(i) + dqMax).*(dqMax - x(i));
    end
    hXFFN{11} = @(x)[[-1,0;0,-1]*PosSwingFN(x(1:5) + 0.0013);
                 -JYSwingFN(x(1:5))*x(6:10)];
%     hXFFN{11} = @(x) x(6) + x(7) + x(8)/2 - 0.5; % Negative theta velocity
    
    hSFN = @(x) [[1,0;0,-1]*PosSwingFN(x(1:5) + 0.0013);
                 -JYSwingFN(x(1:5))*x(6:10)];
    
    qErrMinus = @(x) x(1:5) - q0FN(thetaFN(x(1:5)),x(11)); % Error off the target angles before impact
    qErrPlus = @(x)[1,0,0,0,0;
                    0,0,0,1,0;
                    0,0,0,0,1;
                    0,1,0,0,0;
                    0,0,1,0,0]*qErrMinus(x);
    q0Plus = @(x) q0FN(thetaFN(x([1;4;5;2;3])),x(11)); % Target angles post-impact
    qPlus = @(x) q0Plus(x) + qErrPlus(x);
    
    dqErrMinus = @(x) x(6:10) - dq0FN(thetaFN(x(1:5)),x(11),dthetaFN(x(1:5),x(6:10)),x(12));
    dqErrPlus = @(x) DeltaDqFN(x(1:5))*dqErrMinus(x);
    dq0Plus = @(x) dq0FN(thetaFN(qPlus(x)),x(11),dthetaFN(qPlus(x),DeltaDqFN(x(1:5))*x(6:10)),x(12));
    dqPlus = @(x) dq0Plus(x) + dqErrPlus(x);
    
    RFN  = @(x) [qPlus(x);
                 dqPlus(x);
                 x(11);
                 x(12);
                 x(13)];

    RabbitModel = HybridDynamicModel({nStates},{fFN},{guFN},{[]},{hXFN},{hXFFN},{hSFN},{RFN});

end