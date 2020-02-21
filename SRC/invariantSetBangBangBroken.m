% Could add new su for each hXF

function [out] = invariantSetBangBang(x,f,gu,gd,hX,hXF,vu_0,c_0,uSel,dl,degree,x_TP,TP_weight)
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
    prob.uSel = uSel;

    %% Initialize variables
    v = cell(length(x),1);
    obj = -inf;
    maxIter = 500;
    
    convFlag = 0;

    %% Get initial values for lambda, c, su
    [lambda, c, su] = lOPT(prob, vu_0, 1);
    fprintf('Lambda and u found found for c = %4.4f\n',c)
    c = c_0;
    
    %% Alternate between lambda and v optimization
    for i = 1:maxIter
        % Uncomment this for weighting the integral of v in areas where
        % lambda > 0. This causes trouble with convergence checking
%         [imonom,icoeff] = indFun(lambda,x,hX,dl,degree,options);
        
        c = ~convFlag*max(c,0);

        imonom = msspoly(1); icoeff = 1;

        [vNew, vuNew, qNew, objNew, area] = vOPT(prob, lambda, imonom, icoeff, c, su, i == 1);

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
            vu = vuNew;
            q = qNew;
            if c <= 0
                obj = objNew;
                convFlag = 1;
            end

            [lambda, c, su] = lOPT(prob, vu, 0, q);
            fprintf('Lambda and u found found for c = %4.4f\n',c)

            if vMax < eps
                c = c_0;
            end
        end 
    end

    out.v = v;
    out.vu = vu;
    out.su = su;
    out.q = q;
    out.obj = obj;
    out.lambda = lambda;
    out.c = c;
end


function [vOut,vuOut,qOut,objOut,areaOut] = vOPT(prob, lambda, imonom, icoeff, c, su, init)
    %% v optimization
    %% Unpack problem
    x = prob.x;
    f = prob.f;
    gu = prob.gu;
    uSel = prob.uSel;
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
        
    nd_ = size(gd,2);
    nvu_ = size(uSel,2);
    
    %% Define the persistents
    persistent progSaved_ v_ vu_ vmonom_ vcoeff_ Lfvu_ Lgdvu_ q_ e_
    
    if init
        %% define the program
        progSaved_ = spotsosprog;
        
        %% Define the polynomials
        progSaved_ = progSaved_.withIndeterminate( x );

        vmonom_ = monomials( x, 0:degree );
        [ progSaved_, v_, vcoeff_ ] = progSaved_.newFreePoly( vmonom_ );
        
        vu_ = msspoly(zeros(nvu_,1));
        Lfvu_ = msspoly(zeros(nvu_,1));
        
        if isempty(gd)
            withDist = 0;
            q_ = zeros(1,nvu_);
            Lgdvu_ = q_;
            e_ = q_;
        else
            withDist = 1;
            q_ = msspoly(zeros(nvu_,nd_));
            Lgdvu_ = q_;
            e_ = ones(nd_,1);
            qmonom_ = monomials( x, 0:degree );
        end
        
        for i = 1:nvu_
            [ progSaved_, vu_(i) ] = progSaved_.newFreePoly(vmonom_);
            dvudx_ = diff( vu_(i), x );
            Lfvu_(i) = dvudx_ * f + dvudx_ * gu * uSel(:,i);
            
            if withDist
                for j = 1:nd_
                    [ progSaved_, q_(i,j) ] = progSaved_.newFreePoly( qmonom_ );
                end
                Lgdvu_(i,:) = dvudx_ * gd;
            end
        end

        %% Set the SOS constraints that are consistent across calls

        % v < 1
        progSaved_ = sosOnK(progSaved_, 1 - v_, x, hX, degree);

        % v < vu_
        for i = 1:nvu_
            progSaved_ = sosOnK(progSaved_, vu_(i) - v_, x, hX, degree);
        end

        % Disturbances are outside [-1,1]
        for i = 1:nvu_
            for j = 1:nd_
                progSaved_ = sosOnK(progSaved_, q_(i,j) - Lgdvu_(i,j), x, hX, degree);
                progSaved_ = sosOnK(progSaved_, q_(i,j) + Lgdvu_(i,j), x, hX, degree);
            end
        end
    end
    
    prog = progSaved_;
    
    %% Add the SoS constraint specific to this call
    for i = 1:nvu_
        % dvu_dt > 0 on the barrier {x | vu_i(x) = 0}
        prog = sosOnK(prog,Lfvu_(i) - q_(i,:)*e_ + lambda(i)*vu_(i) + c - epsBar, x, hX, degree);
        
        % min(vu_) < 0 on XF
        for j = 1:length(hXF)
            prog = sosOnK(prog, -sum(su.'*vu_) + c - epsFail, x, [hXF{j};hX], degree);
        end
    end
    
    
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
    vuOut = sol.eval(vu_);
    qOut = sol.eval(q_);
    areaOut = double(sol.eval(area_));
    objOut = -double(sol.eval(obj_));
end

function [lambdaOut, cOut, suOut] = lOPT(prob, vu, init, q)
    %% lambda optimization
    %% Unpack problem
    x = prob.x;
    f = prob.f;
    gu = prob.gu;
    gd = prob.gd;
    hX = prob.hX;
    hXF = prob.hXF;
    uSel = prob.uSel;
    degree = prob.degree;
    options = prob.options;
    epsBar = prob.epsBar;
    epsFail = prob.epsFail;
    
    nvu = size(uSel,2);
    monom = monomials( x, 0:degree );
    
    %% Define the persistents (if undefined)
    persistent progSaved_ e_ lambda_ c_ su_ obj_
    
    if init
        %% define the program
        progSaved_ = spotsosprog;
        
        %% Define the noise scaling constant
        [ progSaved_, c_ ] = progSaved_.newFree( 1, 1 );
%         [ progSaved_, c2_ ] = progSaved_.newFree( nvu, 1 );
        
        %% Add the indeterminate
        progSaved_ = progSaved_.withIndeterminate( x );

        %% Define the polynomials
        lambda_ = msspoly(zeros(nvu,1));
        su_ = msspoly(zeros(nvu,1));
        
        for i_ = 1:nvu
            [ progSaved_, lambda_(i_)] = progSaved_.newFreePoly( monom );
            [ progSaved_, su_(i_)] = progSaved_.newFreePoly( monom );
        end
        
        if isempty(gd)
            e_ = 0;
        else
            e_ = ones( size( gd, 2 ), 1 );
        end
        
        %% Add the SoS constraint specific to this call
        for i_ = 1:nvu
            progSaved_ = sosOnK(progSaved_, su_(i_) , x, hX, degree);
        end
        
        %% Set the objective
        obj_ = c_;
    end
    
    prog = progSaved_;
    
    %% Add q's if not provided
    if nargin < 4
        if isempty(gd)
            q = msspoly(zeros(nvu,1));
        else
            q = msspoly(zeros(nvu,size(gd,2)));
            for i = 1:nvu
                dvudx = diff( vu(i), x );
                Lgdvu = dvudx * gd;
                
                for j = 1:size(gd,2)
                    [prog, q(i,j)] = prog.newFreePoly( monom );
                
                    prog = sosOnK(prog, q(i,j) - Lgdvu(j), x, hX, degree);
                    prog = sosOnK(prog, q(i,j) + Lgdvu(j), x, hX, degree);
                end
            end
        end
    end
    
    
    %% Set the SOS constraints specific to this function call
    for i = 1:nvu
        dvudx = diff( vu(i), x );
        Lfvu = dvudx * f + dvudx * gu * uSel(:,i);
        
        prog = sosOnK(prog,Lfvu - q(i,:)*e_  + lambda_(i)*vu(i) + c_ - epsBar, x, hX, degree);
        
        % min(vu_) < 0 on XF
        for j = 1:length(hXF)
            prog = sosOnK(prog, -sum(su_.'*vu) + c_ - epsFail, x, [hXF{j};hX], degree);
        end
    end

    %% Solve
    sol = prog.minimize( obj_, @spot_mosek, options );

    lambdaOut = sol.eval(lambda_);
    suOut = sol.eval(su_);
    cOut = double(sol.eval(c_));
end