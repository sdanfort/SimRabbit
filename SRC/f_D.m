function f_Dubinesque = f_D(in1,in2)
%F_D
%    F_DUBINESQUE = F_D(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 7.2.
%    06-Apr-2018 17:03:15

uD1 = in2(1,:);
uD2 = in2(2,:);
xS2 = in1(2,:);
xS3 = in1(3,:);
xS4 = in1(4,:);
xS5 = in1(5,:);
t2 = 1.0./xS5;
f_Dubinesque = [xS2+xS3.*xS5;xS3.*1.829239766081871e1-xS4.*2.830409356725146e2-t2.*xS2.*1.39625730994152e2+t2.*xS4.*1.746226900584796e2+uD1.*xS5.*2.830409356725146e2+xS4.*(t2.*2.085333333333333e1-xS5);xS4;xS4.*-2.0e2+uD1.*xS5.*2.0e2;uD2.*6.666666666666667e-4-xS5.^2.*2.003333333333333e-4-1.4715e-1];
