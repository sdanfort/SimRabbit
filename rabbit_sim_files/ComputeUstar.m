function [uStar, q0, dq0, u_0, u_ua, MinvT] = ComputeUstar(xHat,ua)
% Note autogen needs to be added to path before this can work

tauSel = [0,0,0,0;
          1,0,0,0;
          0,1,0,0;
          0,0,1,0;
          0,0,0,1];

if nargin < 1
    ua = 0;
end

theta = xHat(1);
dtheta = xHat(2);
alpha = xHat(3);
dalpha = xHat(4);

q0 = q0FN(theta,alpha);
dq0 = dq0FN(theta,alpha,dtheta,dalpha);

M = MassMatrixFN(q0);
Fcg = CoriGravFN(q0,dq0);
dhdq_ = dhdqFN(q0,alpha);

MinvT = M\tauSel;
MTterm = dhdq_*MinvT;
MFterm = dhdq_*(M\Fcg);

u_0 = -MTterm\(dhLumpFN(q0,alpha,dq0,dalpha) + MFterm);
u_ua = -MTterm\dhdaFN(q0,alpha);
uStar = u_0 + u_ua*ua;

end