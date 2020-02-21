function [out] = invariantSet(x,f,gu,gd,hX,hXF,lambda_0,c_0,u_0,dl,degree,x_TP,TP_weight,du_Const)
%
%  x      -- n-by-1 free msspoly
%  f      -- n-by-1 msspoly in x
%  g      -- n-by-q msspoly in x
%  hX     -- m-by-1 msspoly in x
%  hXF{j} -- p{j}-by-1 msspoly in x
%  dl     -- function mapping n-by-1 msspolys into n-by-1 double.
%  degree -- positive scalar integer
%
%  Computes an inner approximation to the largest invariant set for:
%
%  xdot = f(x) + g(x)*d,
%
%  for all d in D, that doesn't contain XF
%
%  D  = { d | abs(d(i)) < 1}
%
%  XF = UNION{XF_j}, where XF_j = {x | hXF{j}(i) > 0 for all i} 
%
%  dlambda(p) = int_X p(x) dx component-wise (lebesgue moments)
%
%  degree dictates the order of certain approximations.
%
%  r(x) >= 1 is an approximation of x(0) that can stay inside r(x) >=1
%          indefinitely.

    tic

    %% Load input arguments into struct
    prob.x = x;
    prob.f = f;
    prob.gu = gu;
    prob.gd = gd;
    prob.hX = hX;
    prob.hXF = hXF;
    if nargin < 12
        prob.x_TP = [];
        prob.TP_weight = 0;
%         prob.x_TP = double(0*x);
%         prob.TP_weight = 1000;
    else
        prob.x_TP = x_TP;
        prob.TP_weight = TP_weight; 
    end
    if nargin < 14
        prob.du_Const = [];
    else
        prob.du_Const = du_Const;
    end
    prob.dl = dl;
    prob.degree = degree;
    
    options = spot_sdp_default_options();
    options.verbose = 0;
    options.domain_size = 1;
    options.solveroptions = [];
    
    prob.options = options;
    prob.epsBar = 1e-3;
    prob.epsFail = 1e-3;
    prob.nModes = length(x);

    %% Initialize variables
    v = cell(length(x),1);
    lambda = lambda_0;
    
    c = c_0;
    u = u_0;
    obj = -inf;
    maxIter = 500;
    
    convFlag = 0;

    %% Alternate between lambda and v optimization
    for i = 1:maxIter
        % Uncomment this for weighting the integral of v in areas where
        % lambda > 0. This causes trouble with convergence checking
%         [imonom,icoeff] = indFun(lambda,x,hX,dl,degree,options);
        
        c = ~convFlag*max(c,0);

        imonom = msspoly(1); icoeff = 1;

        [vNew, qNew, objNew, area] = vOPT(prob, lambda, imonom, icoeff, c, u, i == 1);

        vMax = polyMax(vNew,x,hX,degree,options);
        
        fprintf('V found for c = %4.2f with objective %4.2f, area %4.2f and max value %4.2f\n',max(c,0),objNew,area,vMax)

        if (c <= 0) && vMax > eps && (objNew < (1+1e-5)*obj) 
            fprintf('Feasible solution found\n')
            fprintf('Elapsed time: %4.2f\n',toc)
%             v = vNew;
%             q = qNew;
%             obj = objNew;
            break;
        else
            v = vNew;
            q = qNew;
            if c <= 0
                obj = objNew;
                convFlag = 1;
            end

            [lambda, c, u] = lOPT(prob, v, q, i == 1);
            fprintf('Lambda and u found found for c = %4.4f\n',c)

            if vMax < eps
                c = c_0;
            end
        end 
    end

    out.v = v;
    out.q = q;
    out.obj = obj;
    out.lambda = lambda;
    out.c = c;
    out.u = u;
end


function [vOut,qOut,objOut,areaOut] = vOPT(prob, lambda, imonom, icoeff, c, u, init)
    %% v optimization
    %% Unpack problem
    x = prob.x;
    f = prob.f;
    gu = prob.gu;
    gd = prob.gd;
    hX = prob.hX;
    hXF = prob.hXF;
    dl = prob.dl;
    x_TP = prob.x_TP;
    TP_weight = prob.TP_weight;
    degree = prob.degree;
    options = prob.options;
    epsBar = prob.epsBar;
    epsFail = prob.epsFail;
    
    %% Define the persistents
    persistent progSaved_ v_ vmonom_ vcoeff_ Lfv_ Lguv_ Lgdv_ q_ e_
    
    if init
        %% define the program
        progSaved_ = spotsosprog;
        
        %% Define the polynomials
        progSaved_ = progSaved_.withIndeterminate( x );

        vmonom_ = monomials( x, 0:degree );
        [ progSaved_, v_, vcoeff_ ] = progSaved_.newFreePoly( vmonom_ );

        % creating the variables that will be used later
        dvdx_ = diff( v_, x );
        Lfv_ = dvdx_ * f;

        if isempty(gu)
            Lguv_ = 0;
        else
            Lguv_ = dvdx_ * gu;
        end

        if isempty(gd)
            q_ = 0;
            Lgdv_ = 0;
            e_ = 0;
        else
            e_ = ones( size( gd, 2 ), 1 );
        
            qmonom_ = monomials( x, 0:degree );
            q_ = msspoly( zeros( [size( gd, 2 ), 1] ) );
            for j = 1:size(gd,2)
                [ progSaved_, q_(j) ] = progSaved_.newFreePoly( qmonom_ );
            end

            Lgdv_ = dvdx_ * gd;
        end

        %% Set the SOS constraints that are consistent across calls

        % v < 1
        progSaved_ = sosOnK(progSaved_, 1 - v_, x, hX, degree);

        % v < 0 on failure set
        for i = 1:length(hXF)
            progSaved_ = sosOnK(progSaved_, 1000*(-v_ - epsFail), x, hXF{i}, degree);
        end

        % v < 0 on boundary of X
        for i = 1:size(hX)
            progSaved_ = sosOnK(progSaved_, 1000*(-v_ - epsFail), x, [hX; -hX(i)], degree);
        end

        % Disturbances are outside [-1,1]
        for i = 1:size( gd, 2 )
            progSaved_ = sosOnK(progSaved_, q_(i) - Lgdv_(i), x, hX, degree);
            progSaved_ = sosOnK(progSaved_, q_(i) + Lgdv_(i), x, hX, degree);
        end
    end
    
    prog = progSaved_;
    
    %% Add the SoS constraint specific to this call
    % dv_idt > 0 on the barrier {x | v_i(x) = 0}
    prog = sosOnK(prog,Lfv_ - e_'*q_ + Lguv_*u + lambda*v_ + c - epsBar, x, hX, degree);
    
    %% Set the objective
    area_ = icoeff' * dl( imonom * vmonom_' ) * vcoeff_;
    
    obj_ = -area_;
    
    if ~isempty(x_TP)
        obj_ = obj_ - TP_weight*sum(msubs(v_,x,x_TP));
    end
        % Possible improvements to the single tentpole idea:
        %   1) Multiple tentpoles (e.g. around the limit cycle)
        %   2) Tentpole with area (e.g. circle or box with integral)
        %   3) Adaptive tentpoles
        %       a) pick tentpoles inside previous v using grid
        %       b) pick largest tent platform that fits in previous v
        %   4) Use tentpoles as a constraint (v(x_TP) > 0) once it goes up
    %% Solve
    sol = prog.minimize( obj_, @spot_mosek, options );

    vOut = sol.eval(v_);
    qOut = sol.eval(q_);
    areaOut = double(sol.eval(area_));
    objOut = -double(sol.eval(obj_));
end

function [lambdaOut, cOut, uOut] = lOPT(prob, v, q, init)
    %% lambda optimization
    %% Unpack problem
    x = prob.x;
    f = prob.f;
    gu = prob.gu;
    gd = prob.gd;
    hX = prob.hX;
    degree = prob.degree;
    options = prob.options;
    epsBar = prob.epsBar;
    
    %% Define the persistents (if undefined)
    persistent progSaved_ e_ lambda_ c_ u_ obj_
    
    if init
        %% define the program
        progSaved_ = spotsosprog;
        
        %% Define the noise scaling constant
        [ progSaved_, c_ ] = progSaved_.newFree( 1 );
        
        %% Add the indeterminate
        progSaved_ = progSaved_.withIndeterminate( x );

        %% Define the polynomials
        lmonom_ = monomials( x, 0:degree );
        [ progSaved_, lambda_ ] = progSaved_.newFreePoly( lmonom_ );

        if isempty(gu)
            u_ = 0;
        else
            pmonom_ = monomials( x, 0:degree);
            u_ = msspoly( zeros( [size( gu, 2 ), 1] ) );
            for j = 1:size(gu,2)
                [ progSaved_, u_(j) ] = progSaved_.newFreePoly( pmonom_ );
            end
        end

        if isempty(gd)
            e_ = 0;
        else
            e_ = ones( size( gd, 2 ), 1 );
        end

        %% Set the SOS constraints consistent across calls

        % Control inputs are inside [-1,1]
        for i = 1:size( gu, 2 )
            progSaved_ = sosOnK(progSaved_, 1 - u_(i), x, hX, degree);
            progSaved_ = sosOnK(progSaved_, 1 + u_(i), x, hX, degree);
        end
        
        %% Set the objective
        obj_ = c_;
    end
    
    prog = progSaved_;
    
    dvdx = diff( v, x );
    Lfv = dvdx * f;

    if isempty(gu)
        Lguv = 0;
    else
        Lguv = dvdx * gu;
    end

    %% Set the SOS constraints specific to this function call
    prog = sosOnK(prog,Lfv - e_'*q + Lguv*u_  + lambda_*v + c_ - epsBar, x, hX, degree);

    %% Solve
    sol = prog.minimize( obj_, @spot_mosek, options );

    lambdaOut = sol.eval(lambda_);
    uOut = sol.eval(u_);
    
    cOut = double(sol.eval(c_));
end