function [out] = invariantSetHybridBang(x,f,gu,gd,sX,R,hX,hXF,lambda_0,sigma_0,c_0,u_0,dl,degree,x_TP,TP_weight,x_Con,plotFun,sExp)
%
%  x{i}         -- n_i-by-1 free msspoly ( \forall i \in {1,...,nModes} )
%  f{i}         -- n_i-by-1 msspoly in x{i} (Defines the dynamcis)
%  gu{i}        -- n_i-by-nu_i msspoly in x{i} (Defines the input)
%  gd{i}        -- n_i-by-nd_i msspoly in x{i} (Defines the disturbance)
%  sX{i,j}      -- l_i-by-1 msspoly in x{i} (defines the guard from i to j)
%  R{i,j}       -- n_j-by-1 msspoly in x{i} (defines the reset from i to j)
%  hX{i}        -- m_i-by-1 msspoly in x{i} (Defines X)
%  hXF{j}{i}    -- p_ji-by-1 msspoly in x{i} (Defines unsafe sets within X)
%  lambda_0{i}  -- 1-by-1 msspoly in x{i} (Initial guess of the invarint set)
%  sigma_0{i,j} -- 1-by-1 msspoly in x{i} (Initial guess of s variable on
%                  the guard from i to j)
%  c_0          -- Scalar constant (Initial guess of the disturbance scaling)
%  u_0{i}       -- nu_i-by-1 msspoly in x{i} (Initial guess of the controller)
%  dl{i}        -- function mapping n_i-by-1 msspoly's into n-by-1 double.
%  degree       -- positive scalar integer
%
%  Computes an inner approximation to the largest invariant that doesn't 
%  contain XF for:
%
%  xdot{i} = f{i}(x) + gd{i}(x)*d{i} + gu{i}(x)*u{i},
%
%  for some u{i} in U{i}, and for all d{i} in D{i}
%
%  D{i}  = { d{i} | abs(d{i}(j)) < 1}
%  U{i}  = { u{i} | abs(d{i}(j)) < 1}
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
    prob.sX = sX;
    prob.R = R;
    prob.hX = hX;
    prob.hXF = hXF;
    if nargin < 15
        prob.x_TP = repcell({[]},length(x),1);
        prob.TP_weight = 0;
    else
        prob.x_TP = x_TP;
        prob.TP_weight = TP_weight; 
    end
    if nargin < 17
        prob.x_Con = repcell({[]},length(x),1);
    else
        prob.x_Con = x_Con;
    end
    prob.dl = dl;
    prob.degree = degree;
    if ~exist('sExp')
        prob.sExp = repcell({[]},nModes,nModes);
    else
        prob.sExp = sExp;
    end
    
    options = spot_sdp_default_options();
    options.verbose = 0;
    options.domain_size = 1;
    options.solveroptions = [];
    
    prob.options = options;
    prob.epsBar = 0*1e-2;
    prob.epsFail = 1e-2;
    prob.nModes = length(x);

    %% Initialize variables
    v = cell(length(x),1);
    lambdaNew = lambda_0;
    sigma = sigma_0;
    
    c = c_0;
    uNew = u_0;
    obj = -inf;
    maxIter = 1000;
    epsConv = 1e-3;
    convFlag = 0;

    %% Alternate between lambda and v optimization
    for i = 1:maxIter
        
%         [imonom,icoeff] = indFun(lambda,x,hX,dl,degree,options);
%         imonom  = repcell({msspoly(1)},nModes,1);
%         icoeff = repcell({1},nModes,1);

        c = ~convFlag*max(c,0);

        [vNew, qNew, rhoNew, objNew, area] = vOPT(prob, lambdaNew, sigma, c, uNew, i == 1);

        vMaxNew = polyMax(vNew,x,hX,degree,options);
        
        fprintf('V found for c = %4.2f with objective %4.2f, area %4.2f and max value %4.4f\n',c,objNew,area,vMaxNew)
        
        toc

        if nargin > 17 && ~isempty(plotFun)
            plotFun(x,vNew,qNew,lambdaNew,uNew)
        end
        
        if (i > 1) && ( (convFlag || (c <= 0)) && (objNew < obj + epsConv*abs(obj)) ) % || (vMaxNew <= epsConv) ) )
            fprintf('Feasible solution found\n')
%             v = vNew;
%             q = qNew;
%             obj = objNew;
%             u = uNew;
%             lambda = lambdaNew;
            break;
        else            
            v = vNew;
            q = qNew;
            u = uNew;
            rho = rhoNew;
            lambda = lambdaNew;
            if c <= 0
                obj = objNew;
                vMax = vMaxNew;
                if vMax > eps
                    convFlag = 1;
                end
            end
            
            [lambdaNew, cNew, uNew] = lOPT(prob, v, q, rho, i == 1);
            fprintf('Lambda sigma and u found found for c = %4.4f\n',cNew)

            toc
            
            if vMaxNew < eps
                c = c_0;
            else
                c = min(c,cNew);
            end
        end 
    end

    out.v = v;
    out.q = q;
    out.rho = rho;
    out.obj = obj;
    out.lambda = lambda;
    out.sigma = sigma;
    out.c = c;
    out.u = u;
end


function [vOut,qOut,rhoOut,objOut,areaOut] = vOPT(prob, lambda, sigma, c, u, init)
    %% v optimization
    
    %% Unpack problem
    x = prob.x;
    f = prob.f;
    gu = prob.gu;
    gd = prob.gd;
    sX = prob.sX;
    sExp = prob.sExp;
    R = prob.R;
    hX = prob.hX;
    hXF = prob.hXF;
    dl = prob.dl;
    x_TP = prob.x_TP;
    TP_weight = prob.TP_weight;
    x_Con = prob.x_Con;
    degree = prob.degree;
    options = prob.options;
    epsFail = prob.epsFail;
    epsBar = prob.epsBar;
    nModes = prob.nModes;
    
    %% Define the persistents
    persistent progSaved v Lfv Lguv Lgdv q e obj TP_val area uSel rho
    
    if init
        %% define the program
        progSaved = spotsosprog;

        TP_val = 0;
        area = 0;
        v = cell(nModes,1);
        Lfv = cell(nModes,1);
        Lguv = cell(nModes,1);
        Lgdv = cell(nModes,1);
        q = cell(nModes,1);
        e= cell(nModes,1);
        uSel = cell(nModes,1);
        rho = cell(nModes,1);
        
        %% Define the polynomials
        for i = 1:nModes
            progSaved = progSaved.withIndeterminate( x{i} );
            
            monom_ = monomials( x{i}, 0:degree );
            [ progSaved, v{i}, vcoeff_ ] = progSaved.newFreePoly( monom_ );

            % creating the variables that will be used later
            dvdx = diff( v{i}, x{i} );
            Lfv{i} = dvdx * f{i};

            if isempty(gu{i})
                Lguv{i} = 0;
                u{i} = 0;
                uSel{i} = 0;
                rho{i} = 0;
            else
                nu_ = size(gu{i},2);
                Lguv{i} = dvdx * gu{i};
                uSel{i} = nchoosek(repmat([-1,1],1,nu_),nu_).';
                rho{i} = msspoly(zeros(nu_,size(uSel,2)));
                for j = 1:nu_
                    for k = 1:size(uSel{i},2)
                        [ progSaved, rho{i}(j,k) ] = progSaved.newFreePoly( monom_ );
                    end
                end
            end

            if isempty(gd{i})
                q{i} = 0;
                Lgdv{i} = 0;
                e{i} = 0;
            else
                e{i} = ones( size( gd{i}, 2 ), 1 );
                
                q{i} = msspoly( zeros( [size( gd{i}, 2 ), 1] ) );
                for j = 1:size(gd{i},2)
                    [ progSaved, q{i}(j) ] = progSaved.newFreePoly( monom_ );
                end

                Lgdv{i} = dvdx * gd{i};
            end

            %% Set the SOS constraints that are consistent across calls

            % v < 1
            progSaved = sosOnK(progSaved, 1 - v{i}, x{i}, hX{i}, degree);

            % v < 0 on failure set
            for j = 1:length(hXF{i})
                progSaved = sosOnK(progSaved, 100*(-v{i} - epsFail), x{i}, [hX{i};hXF{i}{j}], degree);
            end
            
            for j = 1:size( gu{i}, 2 )
                for k = 1:size(rho{i},2)
                    progSaved = sosOnK(progSaved, rho{i}(j,k), x{i}, hX{i}, degree);
                end
            end
                    
            % Disturbances are outside [-1,1]
            for j = 1:size( gd{i}, 2 )
                progSaved = sosOnK(progSaved, q{i}(j) - Lgdv{i}(j), x{i}, hX{i}, degree);
                progSaved = sosOnK(progSaved, q{i}(j) + Lgdv{i}(j), x{i}, hX{i}, degree);
            end

            %% Set the objective
            area = area + dl{i}( monom_' ) * vcoeff_;
            if ~isempty(x_TP{i})
                TP_val = TP_val + TP_weight*sum(msubs(v{i},x{i},x_TP{i}));
            end
            
            if ~isempty(x_Con{i})
                for j = 1:size(x_Con{i},2)
                    progSaved = withPos(progSaved, msubs(v{i},x{i},x_Con{i}(:,j)) - epsFail);
                end
            end
        end
        obj = -TP_val - area;
    end
    
    prog = progSaved;
    
    %% Add the SoS constraint specific to this call
    for i = 1:nModes
        % dv_idt > 0 on the barrier {x | v_i(x) = 0}
        for j = 1:size(uSel{i},2)
            prog = sosOnK(prog,Lfv{i} - e{i}'*q{i} + Lguv{i}*uSel{i}(:,j) + ...
                               lambda{i}*v{i} - rho{i}(:,j).'*(uSel{i}(:,j).*u{i}) + ...
                               c - epsBar, x{i}, [hX{i}], degree);
        end

        for j = 1:nModes
            if ~isempty(sExp{j,i})
                hX_ = subs(hX{j},x{j},sExp{j,i});
                R_ = subs(R{j,i},x{j},sExp{j,i});
                sigma_ = subs(sigma{j,i},x{j},sExp{j,i});
                vIn_ = subs(v{j},x{j},sExp{j,i});
                vOut_ = subs(v{i},x{i},R_);
                x_ = decomp(sExp{j,i});
                
                prog = sosOnK(prog, 100*(vOut_ - sigma_*vIn_), x_, hX_, degree);
                
            elseif ~isempty(sX{j,i}) % If there is a guard from mode j to mode i
                % v_i(R_ji(x)) >= 0 on the set {x | v_j(x) >= 0, S_j(x) = 0)
                prog = sosOnK(prog, 100*(subs(v{i},x{i},R{j,i}) - sigma{j,i}*v{j} - epsBar), x{j}, [hX{j}; sX{j,i}], degree);
            end
        end
    end
    
    %% Solve
    sol = prog.minimize( obj, @spot_mosek, options );

    vOut = cell(nModes,1);
    qOut = cell(nModes,1);
    rhoOut = cell(nModes,1);
    for i = 1:nModes
        vOut{i} = sol.eval(v{i});
        qOut{i} = sol.eval(q{i});
        rhoOut{i} = sol.eval(rho{i});
    end
    areaOut = double(sol.eval(area));
    objOut = -double(sol.eval(obj));
end

function [lambdaOut, cOut, uOut] = lOPT(prob, v, q, rho, init)
    %% lambda optimization
    %% Unpack problem
    x = prob.x;
    f = prob.f;
    gu = prob.gu;
    gd = prob.gd;
    hX = prob.hX;
    degree = prob.degree;
    options = prob.options;
    nModes = prob.nModes;
    epsBar = prob.epsBar;
    
    %% Define the persistents (if undefined)
    persistent progSaved e lambda c u obj uSel
    
    if init
        e = cell(nModes,1);
        lambda = cell(nModes,1);
        u = cell(nModes,1);
        
        %% define the program
        progSaved = spotsosprog;
        
        %% Define the noise scaling constant
        [ progSaved, c ] = progSaved.newFree( 1 );
%         c = 0;

        for i = 1:nModes
            %% Add the indeterminate
            progSaved = progSaved.withIndeterminate( x{i} );

            %% Define the polynomials
            lmonom = monomials( x{i}, 0:degree );
            [ progSaved, lambda{i} ] = progSaved.newFreePoly( lmonom );

            if isempty(gu{i})
                u{i} = 0;
            else
                nu_ = size(gu{i},2);
                lmonom = monomials( x{i}, 0:degree);
                u{i} = msspoly( zeros( [nu_, 1] ) );
                for j = 1:nu_
                    [ progSaved, u{i}(j) ] = progSaved.newFreePoly( lmonom );
                end
                uSel{i} = nchoosek(repmat([-1,1],1,nu_),nu_).';
            end

            if isempty(gd{i})
                e{i} = 0;
            else
                e{i} = ones( size( gd{i}, 2 ), 1 );
            end

            %% Set the SOS constraints consistent across calls

            % Control inputs are inside [-1,1]
%             for j = 1:size( gu{i}, 2 )
%                 progSaved = sosOnK(progSaved, 1 - u{i}(j), x{i}, hX{i}, degree);
%                 progSaved = sosOnK(progSaved, 1 + u{i}(j), x{i}, hX{i}, degree);
%             end
        end
        
        %% Set the objective
        obj = c;
    end
    
    prog = progSaved;
    
    for i = 1:nModes
        dvdx = diff( v{i}, x{i} );
        Lfv = dvdx * f{i};

        if isempty(gu{i})
            Lguv = 0;
        else
            Lguv = dvdx * gu{i};
        end

        %% Set the SOS constraints specific to this function call
        for j = 1:size(uSel{i},2)
            prog = sosOnK(prog,Lfv - e{i}'*q{i} + Lguv*uSel{i}(:,j) + ...
                               lambda{i}*v{i} - rho{i}(:,j).'*(uSel{i}(:,j).*u{i}) + ...
                               c - epsBar, x{i}, [hX{i}], degree);
        end
    end
    
%     prog = withPos(prog, cOld - c);

    %% Solve
    sol = prog.minimize( obj, @spot_mosek, options );

    lambdaOut = cell(nModes,1);
    uOut = cell(nModes,1);
    
    for i = 1:nModes
        lambdaOut{i} = sol.eval(lambda{i});
        uOut{i} = sol.eval(u{i});
    end
    
    cOut = double(sol.eval(c));
end
