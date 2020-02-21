function [dThetaPlusVals] = sampleReset(xHatMinusVals)
    nPoints = size(xHatMinusVals,1);
    dThetaPlusVals = zeros(nPoints,1);
    
    for i = 1:size(xHatMinusVals,1)
        thetaMinus_ = xHatMinusVals(i,1);
        dthetaMinus_ = xHatMinusVals(i,2);
        alphaMinus_ = xHatMinusVals(i,3);
        dalphaMinus_ = xHatMinusVals(i,4);
        
        qMinus_ = q0FN( thetaMinus_, alphaMinus_ );
        qPlus_ = [ 1,0,0,0,0;
                    0,0,0,1,0;
                    0,0,0,0,1;
                    0,1,0,0,0;
                    0,0,1,0,0 ]*qMinus_;
        dqMinus_ = dq0FN( thetaMinus_, alphaMinus_, dthetaMinus_,...
            dalphaMinus_ );
        dqPlus_ = DeltaDqFN(qMinus_)*dqMinus_;
        dThetaPlusVals(i) = dthetaFN(qPlus_,dqPlus_);
    end
end