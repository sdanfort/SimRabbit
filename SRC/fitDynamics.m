function [f,gu,gd] = fitDynamics( xVal, uVal, dxVal, dxNomVal, fFitBasis, guFitBasis, gdMax, w, fitInds, x, hX, dl, degree )
%BOUNDPOINTS Find the smallest positive polynomial that bounds the absolute
%            value of the input points

f = msspoly(zeros(length(x),1));
gu = msspoly(zeros(size(guFitBasis)));
gd = msspoly(zeros(length(x),length(fitInds)));
nu = size(guFitBasis,2);

for i = fitInds
    prog = spotsosprog;

    options = spot_sdp_default_options();
    options.verbose = 1;
    options.domain_size = 1;
    options.solveroptions = [];

    prog = prog.withIndeterminate( x );

    
    % Define dynamics objects
    [ prog, fFit ] = prog.newFreePoly( fFitBasis{i} );
    guFit = msspoly(zeros(1,nu));
    for j = 1:nu
        [ prog, guFit(j) ] = prog.newFreePoly( guFitBasis{i,j} );
    end
    gdMonom = monomials( x, 0:degree );
    [ prog, gdFit, gcoeff_ ] = prog.newFreePoly( gdMonom );

    errVal = dxVal(i,:) - ( dxNomVal(i,:) + msubs(fFit,x,xVal) + sum(msubs(guFit',x,xVal).*uVal,1));

    gdVal = msubs(gdFit,x,xVal);
    
    prog = withPos(prog, gdVal - errVal);
    prog = withPos(prog, gdVal + errVal);
    
    prog = sosOnK(prog, gdFit, x, hX, degree);
    prog = sosOnK(prog, gdMax - gdFit, x, hX, degree);
    
    Area = double(dl(msspoly(1)));
    nPoints = size(xVal,2);
    
    c = -w/Area * dl( gdMonom' )*gcoeff_ + 1/nPoints * sum(msubs(gdFit,x,xVal));
%     c = dl(gdMonom')*gcoeff_;
%     c = sum(msubs(gdFit,x,xVal));

    sol = prog.minimize( c, @spot_mosek, options );

    f(i) = sol.eval(fFit);
    gu(i,:) = sol.eval(guFit);
    gd(i,find(fitInds==i,1)) = sol.eval(gdFit);
end

