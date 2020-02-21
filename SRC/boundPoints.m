function g = boundPoints( xVal, yVal, x, hX, dl, degree )
%BOUNDPOINTS Find the smallest positive polynomial that bounds the absolute
%            value of the input points

prog = spotsosprog;

options = spot_sdp_default_options();
options.verbose = 1;
options.domain_size = 1;
options.solveroptions = [];

prog = prog.withIndeterminate( x );

gMonom = monomials( x, 0:degree );
[ prog, g, gcoeff_ ] = prog.newFreePoly( gMonom );

absYVal = abs(yVal);

prog = withPos(prog, msubs(g,x,xVal) - absYVal);

prog = sosOnK(prog, g, x, hX, degree);

c = dl( gMonom' ) * gcoeff_;

sol = prog.minimize( c, @spot_mosek, options );

g = sol.eval(g);

end
