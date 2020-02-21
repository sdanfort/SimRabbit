function [tSim, xSim, uSim, xHatSim, uASim] = SimRabbitAUG(tRange, xHat0,...
    uaFN, plotFlag, qpFlag, Tmax)

% In this function I will reference equations and sections from
% biped_safety_regulation.pdf
% (included in the RabbitPitch folder)

if nargin < 4
    plotFlag = 1;
end
if nargin < 5
    qpFlag = 0;
end
if nargin < 6
    Tmax = 75;
end

%PD controller gains
% kpHZD = 1000;
% kdHZD = 100;

kpHZD = 4000;
kdHZD = 100;

% Define initial manifold parameters from the initial manifold state xHat0
% (Section III-B)
th0 = xHat0(1);
a0 = xHat0(3);
dth0 = xHat0(2);
da0 = xHat0(4);

% Our zero dynamics manifold is given by Z = (q, a, dq, da) (Thm 1)
% q0 evaluated at th0, a0 gives us a 5x1 array
% dq0FN evaluated at th0, a0, dth0, da0 gives us a 5x1 array
% Then we stack a0 and da0 on the bottom
% The last line is time t, which starts at 0.
x0 = [q0FN(th0,a0);
      dq0FN(th0,a0,dth0,da0);
      a0;
      da0;
      0];

% Define the xHat function handle. 
% xHat = [ theta dtheta alpha dalpha ]'
% theta is an autogenerated function of q0, which is x(1:5)
% detheta is an autogenerated function of q0 (x(1:5)) and dq0 (x(6:10))
xHat = @(x) [ thetaFN( x(1:5) );
              dthetaFN( x(1:5), x(6:10) );
              x(11);
              x(12) ];

% Function handle for computing HZD input, this function is then used in
% RabbitModel.Simulate
uFN = @(x) ComputeUHZD(x);

RabbitModel = GenRabbitAug(40);

% Simulate rabbit model from x0.
[ tSim, xAugSim ] = RabbitModel.Simulate( {uFN}, {[]}, tRange, x0, 1 );
% Actual Rabbit states are the first ten states of xAugSim (q0 and dq0)
% So we define these as xSim.
xSim = xAugSim{1}( :, 1:10 );

% Now we need to get our simulated manifold parameter states xHatSim
xHatSim = zeros( length(tSim), 4 );
% And we also want the joint torques.
uAugSim = zeros( length(tSim), 5 );
% uAugSim = zeros( length(tSim), 6 );
for i = 1:length(tSim)
    
    % Evaluating our xHat function at each time step to get the 4 manifold
    % parameters:
    xHatSim(i,:) = xHat(xAugSim{1}(i,:).').';
    
    % Getting the input at each time step with the augmented x-state.
    % ----> Why do we calculate the input AFTER we simulate? Shouldn't the
    % input already be calculated somewhere else?
    uAugSim(i,:) = ComputeUHZD(xAugSim{1}(i,:).').';
    
end

% Now we take only the first four entries to give us our uSim
uSim = uAugSim( :, 1:4 );
% And there is also this uAsim which this function outputs, but is never
% used in TestTargets.
uASim = uAugSim( :, 5 );

% % Now we take only the first four entries to give us our uSim
% uSim = uAugSim( :, 1:5 );
% % And there is also this uAsim which this function outputs, but is never
% % used in TestTargets.
% uASim = uAugSim( :, 6 );

% Plotting
if plotFlag
    figure(1),clf,plot(tSim,hFN(xAugSim{1}(:,1:5).',xAugSim{1}(:,11).'))

    figure(2),clf,plot(thetaFN(xAugSim{1}(:,1:5).'),...
        dthetaFN(xAugSim{1}(:,1:5).',xAugSim{1}(:,6:10).'))

    figure(3),clf,plot(tSim,uAugSim(:,1:4)),...
        legend({'Stance Hip', 'Stance Knee', 'Swing Hip', 'Swing Knee'})

%     animateRabbit(5/max(tSim)*tSim,xAugSim{1},@(x) LinkPosFN(x(1:5)));
end
function [u_] = ComputeUHZD(x_)
    
    % Define our xHat-array given current state x_
    xH_ = xHat(x_);
    
    % Time vector is the last entry of x
    t_ = x_(13);
    
    % Evaluate our shaping parameter alpha using the augmented state,
    % manifold state, and time (see question in TestTargets)
    uaBar_ = uaFN( x_, xH_, t_ );
    
    % We use our current augmented state x_, the state of alpha 
    % (xHat(3) and xHat(4)), and our shaping parameter controller to 
    % generate u_star (the fb controller to stay on the manifold)
    % which is a combination of u_0 and u_ua.
    [ uSt_, u_0_, u_ua_ ] = ComputeUstarFull( x_, [xH_(3); xH_(4)], uaBar_ );
    
    if qpFlag
        
        %Solving the QP in equation 27 for a u_st that respects actuator
        %limits.
        
        H_ = 1;
        f_ = -uaBar_;
        A_ = [u_ua_;-u_ua_];
        b_ = [Tmax - u_0_;
              Tmax + u_0_];
        [ua_,~,exitFlag_] = quadprog( H_, f_, A_, b_, [], [], [], [], uaBar_ );
        
        if exitFlag_ > 0
            uSt_ = ua_;
        end
        
    end
    
    %then... just take take the max between the PD controller and the torque
    %limits
    
    q0_ = q0FN(xH_(1),xH_(3));
    dq0_ = dq0FN(xH_(1),xH_(3),xH_(2),xH_(4));
    
    % The following is our HZD controller plus PD control on 4 of rabbit
    % states.
    u_ = [max(-Tmax,...
          min( Tmax,...
               uSt_ + kpHZD*(q0_(2:5)-x_(2:5)) + kdHZD*(dq0_(2:5)-x_(7:10))));
          uaBar_];

    %ignoring torque limits for now, to make debugging easier:
%     u_ = [ uSt_ + kpHZD*(q0_(2:5)-x_(2:5)) + kdHZD*(dq0_(2:5)-x_(7:10));
%           uaBar_ ];
%      

end

end