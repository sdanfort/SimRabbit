function out1 = YSwingFN(in1)
%YSWINGFN
%    OUT1 = YSWINGFN(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    06-Feb-2019 15:14:22

alphaL = in1(4,:);
alphaR = in1(2,:);
betaL = in1(5,:);
betaR = in1(3,:);
phi = in1(1,:);
out1 = cos(alphaL+phi).*(-2.0./5.0)+cos(alphaR+phi).*(2.0./5.0)-cos(alphaL+betaL+phi).*(2.0./5.0)+cos(alphaR+betaR+phi).*(2.0./5.0);