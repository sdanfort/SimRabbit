function SimRabbitFULL(tRange)
addpath('autogen')

kpHZD = 1000;
kdHZD = 100;
Tmax = 75;

th0 = 0;
a0 = 0.2;
dth0 = 0.5;
da0 = 0;

x0 = [q0FN(th0,a0);
      dq0FN(th0,a0,dth0,da0)];
  
xHat = @(x) [thetaFN(x(1:5));
             a0;
             dthetaFN(x(1:5),x(6:10));
             da0];

uFN = @(x) ComputeUHZD(x);

RabbitModel = GenRabbit(40);

[tSim,xSim] = RabbitModel.Simulate({uFN},{[]},[0,20],x0,1);


figure(1),clf,plot(tSim,hFN(xSim{1}(:,1:5).',a0))

figure(2),clf,plot(thetaFN(xSim{1}(:,1:5).'),dthetaFN(xSim{1}(:,1:5).',xSim{1}(:,6:10).'))

uSim = zeros(length(tSim),4);
for i = 1:length(tSim)
    uSim(i,:) = ComputeUHZD(xSim{1}(i,:).').';
end

figure(3),clf,plot(tSim,uSim),legend({'Stance Hip', 'Stance Knee', 'Swing Hip', 'Swing Knee'})
    
animateRabbit(1*tSim,xSim{1},@(x) LinkPosFN(x(1:5)));

function u_ = ComputeUHZD(x_,ua_)

    [uSt_] = ComputeUstarFull(x_,[a0;0],ua_);
    xH_ = xHat(x_);
    q0_ = q0FN(xH_(1),xH_(2));
    dq0_ = dq0FN(xH_(1),xH_(2),xH_(3),xH_(4));
    
    u_ = max(-Tmax,min(Tmax,uSt_ + kpHZD*(q0_(2:5)-x_(2:5)) + kdHZD*(dq0_(2:5)-x_(7:10))));

end

end