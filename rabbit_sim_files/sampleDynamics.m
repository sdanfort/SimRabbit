function [failVals,FdThetaVals,GdThetaVals,FdAlphaVals,GdAlphaVals] = sampleDynamics(xHatVals,Tmax,Tgap)
    nPoints = size(xHatVals,1);
    FdThetaVals = nan(nPoints,1);
    GdThetaVals = nan(nPoints,1);
    FdAlphaVals = nan(nPoints,1);
    GdAlphaVals = nan(nPoints,1);
    failVals = zeros(nPoints,1);
    
    tauSel = [0,0,0,0;
              1,0,0,0;
              0,1,0,0;
              0,0,1,0;
              0,0,0,1];

    for i = 1:size(xHatVals,1)
        [~, q0_, dq0_, u0_, uua_, ~] = ComputeUstar(xHatVals(i,:),0);
        ua_max = min((sign(uua_)*Tmax - u0_)./uua_);
        ua_min = max((-sign(uua_)*Tmax - u0_)./uua_);
        
        %failVals is positive --> ua does not lie in T limits.
        failVals(i) = Tgap - (ua_max - ua_min);
        
        %if our x_hat does lie within torque limits, though, we sample f and g.
        if  failVals(i) < 0
            M_ = MassMatrixFN(q0_);
            fcg_ = CoriGravFN(q0_,dq0_);
            dtdq_ = dtdqFN(q0_);
            FdThetaVals(i) = dtdq_*(M_\(fcg_ + tauSel*(u0_+uua_*(ua_max+ua_min)/2)));
            GdThetaVals(i) = dtdq_*(M_\(tauSel*uua_*(ua_max-ua_min)/2));
        
            FdAlphaVals(i) = (ua_max+ua_min)/2;
            GdAlphaVals(i) = (ua_max-ua_min)/2;
        end
    end

end