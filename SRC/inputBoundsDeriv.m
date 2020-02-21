function [uPOut, uMOut] = inputBoundsDeriv(x,v,uSafe,f,gu,gd,hX,dl,dudtMax,degree)
    monom = monomials( x, 0:degree );
    
    nu = size( gu, 2 );
    nd = size( gd, 2 );
    
    uPOut = boundOpt(1);
    uMOut = boundOpt(-1);
    
    function uModOut = boundOpt(plusminus)
        options = spot_sdp_default_options();
        options.verbose = 1;
        options.domain_size = 1;
        options.solveroptions = [];
    
        % input plusminus = 1 if positive modification, = -1 if negative
        uModOut = msspoly( zeros(nu, 1) );
        
        for i = 1:nu
            if ~isinf(dudtMax(i))
                % There is a derivative constraint on u(i)
                prog = spotsosprog;
                prog = prog.withIndeterminate( x );

                qu = msspoly( zeros(nu,1) );
                qd = msspoly( zeros(nd,1) );

                [ prog, uMod, uMcoeff ] = prog.newFreePoly( monom );
%                 [prog, sV] = prog.newFreePoly(monom); % S variable for v boundary constraint
                
                u = uSafe + plusminus*uMod;
                
                dudx = diff( u, x );
                Lfu = dudx * f;
                Lguu = dudx * gu;
                if ~isempty(gd)
                    Lgdu = dudx * gd;
                end

                for j = 1:nu
                    [ prog, qu(j) ] = prog.newFreePoly( monom );

                    prog = sosOnK(prog, qu(j) - Lguu(j), x, hX, degree);
                    prog = sosOnK(prog, qu(j) + Lguu(j), x, hX, degree);
                end

                for j = 1:nd
                    [ prog, qd(j) ] = prog.newFreePoly( monom );

                    prog = sosOnK(prog, qd(j) - Lgdu(j), x, hX, degree);
                    prog = sosOnK(prog, qd(j) + Lgdu(j), x, hX, degree);
                end

                prog = sosOnK(prog,  Lfu - sum(qu) - sum(qd) + dudtMax(i), x, [hX;v], degree);
                prog = sosOnK(prog, -Lfu - sum(qu) - sum(qd) + dudtMax(i), x, [hX;v], degree);

                % Mod is less than 0 for x outside safe set
                prog = sosOnK(prog, -uMod, x, [hX;-v], degree);

                obj = - dl( monom' ) * uMcoeff; %-subs(uMod,x,[0;0;0]);

                sol = prog.minimize( obj, @spot_mosek, options );

                uModOut(i) = sol.eval(uMod);
            end
        end
    end
end

