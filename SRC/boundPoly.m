function g = boundPoly( f, x, hX, dl, degree, minCoeff )
%BOUNDERROR Find the smallest polynomial that bounds the absolute
%           value of the input polynomial

prog = spotsosprog;

options = spot_sdp_default_options();
options.verbose = 0;
options.domain_size = 1;
options.solveroptions = [];

prog = prog.withIndeterminate( x );

gMonom = monomials( x, 0:degree );
[ prog, g, gcoeff_ ] = prog.newFreePoly( gMonom );

for i = 1:length(f)
    prog = sosOnK(prog, g - f(i), x, hX, degree);
    prog = sosOnK(prog, g + f(i), x, hX, degree);
end

c = dl( gMonom' ) * gcoeff_;

sol = prog.minimize( c, @spot_mosek, options );

g = sol.eval(g);

if nargin > 5
    [a,b,c] = decomp(g);
    bNew = [];
    cNew = [];
    for i = 1:length(c)
        if abs(c(i)) > minCoeff
            bNew(end+1,:) = b(i,:);
            cNew(1,end+1) = c(i);
        end
    end
    if isempty(cNew)
        g = 0;
    else
        g = recomp(a,bNew,cNew);
    end
end

end

