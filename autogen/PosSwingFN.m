function Swing_Foot_Pos = PosSwingFN(in1)
%POSSWINGFN
%    SWING_FOOT_POS = POSSWINGFN(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    06-Feb-2019 15:14:22

alphaL = in1(4,:);
alphaR = in1(2,:);
betaL = in1(5,:);
betaR = in1(3,:);
phi = in1(1,:);
t2 = alphaL+phi;
t3 = alphaR+phi;
t4 = alphaL+betaL+phi;
t5 = alphaR+betaR+phi;
Swing_Foot_Pos = [sin(t2).*(2.0./5.0)-sin(t3).*(2.0./5.0)+sin(t4).*(2.0./5.0)-sin(t5).*(2.0./5.0);cos(t2).*(-2.0./5.0)+cos(t3).*(2.0./5.0)-cos(t4).*(2.0./5.0)+cos(t5).*(2.0./5.0)];
