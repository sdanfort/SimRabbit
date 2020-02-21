function zOut = findInterpDeriv(x,v,f,gu,gd,hX,dl,dudtMax,degree)
    monom = monomials( x, 0:degree );
    
    nu = size( gu, 2 );
    nd = size( gd, 2 );
    
    zOut = msspoly( zeros(nu, 1) );
    
    options = spot_sdp_default_options();
    options.verbose = 1;
    options.domain_size = 1;
    options.solveroptions = [];
    
    for i = 1:nu
        if ~isinf(dudtMax(i))
            % There is a derivative constraint on u(i)
            prog = spotsosprog;
            prog = prog.withIndeterminate( x );
            
            qu = msspoly( zeros(nu,1) );
            qd = msspoly( zeros(nd,1) );
            
            [ prog, z, zcoeff ] = prog.newFreePoly( monom );
            
            dzdx = diff( z, x );
            Lfz = dzdx * f;
            Lguz = dzdx * gu;
            if ~isempty(gd)
                Lgdz = dzdx * gd;
            end
            
            for j = 1:nu
                [ prog, qu(j) ] = prog.newFreePoly( monom );
                
                prog = sosOnK(prog, qu(j) - Lguz(j), x, hX, degree);
                prog = sosOnK(prog, qu(j) + Lguz(j), x, hX, degree);
            end
            
            for j = 1:nd
                [ prog, qd(j) ] = prog.newFreePoly( monom );

                prog = sosOnK(prog, qd(j) - Lgdz(j), x, hX, degree);
                prog = sosOnK(prog, qd(j) + Lgdz(j), x, hX, degree);
            end
            
            prog = sosOnK(prog,  Lfz - sum(qu) - sum(qd) + dudtMax(i), x, [hX;v], degree);
            prog = sosOnK(prog, -Lfz - sum(qu) - sum(qd) + dudtMax(i), x, [hX;v], degree);
            
            prog = sosOnK(prog, z + 1, x, hX, degree);
            prog = sosOnK(prog, z - 1, x, [hX;-v], degree);
            
            obj = dl( monom' ) * zcoeff;
            
            sol = prog.minimize( obj, @spot_mosek, options );
            
            zOut(i) = sol.eval(z);
        end
    end
end