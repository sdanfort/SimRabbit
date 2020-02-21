function vMax = polyMax(v,x,hX,degree,options)
    if nargin < 5
        options = spot_sdp_default_options();
        options.verbose = 0;
        options.domain_size = 1;
        options.solveroptions = [];
    end

    if ~iscell(v)
        prog = spotsosprog;

        prog = prog.withIndeterminate( x );

        [ prog, c ] = prog.newFree( 1 );

        prog = sosOnK(prog, c - v, x, hX, degree);

        sol = prog.minimize( c, @spot_mosek, options );

        vMax = double(sol.eval(c));
    else
        nModes = length(v);
        
        vMax = -inf;
        
        for i = 1:nModes
            prog = spotsosprog;

            prog = prog.withIndeterminate( x{i} );

            [ prog, c ] = prog.newFree( 1 );

            prog = sosOnK(prog, c - v{i}, x{i}, hX{i}, degree);

            sol = prog.minimize( c, @spot_mosek, options );

            vMax = max(vMax,double(sol.eval(c)));
        end
    end
    
end