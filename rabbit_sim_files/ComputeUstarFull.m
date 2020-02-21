function [uStar,u_0,u_ua] = ComputeUstarFull(x, xa, ua)
% Note autogen needs to be added to path before this can work

% In this function I will reference equations and sections from
% biped_safety_regulation.pdf
% (included in the RabbitPitch folder)
% 
tauSel = [0,0,0,0;
          1,0,0,0;
          0,1,0,0;
          0,0,1,0;
          0,0,0,1];

if nargin < 1
    ua = 0;
end

alpha = xa(1);
dalpha = xa(2);

q = x(1:5);
dq = x(6:10);

M = MassMatrixFN(q);
Fcg = CoriGravFN(q,dq);
dhdq_ = dhdqFN(q,alpha);
% dhLumpFN(q,alpha,dq,dalpha);

MinvT = M\tauSel;
MTterm = dhdq_*MinvT;
MFterm = dhdq_*(M\Fcg);

%These are drawn from paper Eq. 6
u_0 = -MTterm\(dhLumpFN(q,alpha,dq,dalpha) + MFterm);
u_ua = -MTterm\dhdaFN(q,alpha);

uStar = u_0 + u_ua*ua;

end