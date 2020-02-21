function sim_rabbit_gait

kp = 2000;  %this is a huge gain! Can play around with making it smaller 
kd = 100;

%pick an index between 1 and 46. Corresponds to step length values:
%larger index = larger step width
idx = 30;

%load the specific gait from our FROST gait library
[ thSim, xSim, dxSim, x0 ] = load_gait_library( idx );

%fit xSim and dxSim as polynomials of theta
deg = 5;
p = cell(5, 1);
dp = cell(5, 1);
for i = 1:5
    p{i} = polyfit( thSim, xSim(:, i), deg );
    dp{i} = polyfit( thSim, dxSim(:, i), deg );
end

%Now, generate a Rabbit object for simulation
RabbitModel = GenRabbit(40);

%define our input (in function below)
uFN = @(x) u_function( x ); %input function
tf = 5; %final time
m0 = 1; %initial mode
dFN = []; %disturbance (you'll also have to update "GenRabbit" if you add this in)
[ tSim, xSim ] = RabbitModel.Simulate( {uFN}, {dFN}, [0, tf], x0, m0 );

%Animate Rabbit! This function draws the links based on your simulation result
animateRabbit( tSim, xSim{:}, @LinkPosFN );

function u_ = u_function( x_ )
      
    %find our current hip angle:
    th_ = -x_(1) - x_(2) - x_(3)/2;
    
    %now define our target state based on that value:
    q0_ = zeros( 10, 1 );
    for ii = 1:5
        q0_( ii ) = polyval( p{ ii }, th_ );
        q0_( 5 + ii ) = polyval( dp{ ii }, th_ );
    end
    
    %pd control for input:
    u_ = kp*( q0_(2:5)-x_(2:5) ) + kd*( q0_(7:10)-x_(7:10) );

end

end