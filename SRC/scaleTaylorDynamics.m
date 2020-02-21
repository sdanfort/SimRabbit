function [f, gU, gD, scaleParam] = scaleTaylorDynamics(x_US_Sym, x_SC_Poly, f_US, gU_US, gD_US, xRange, uRange, dRange, degree)
% Scale f, gD and gU such that the scaled xRange, uRange and dRange are all
% unit boxes, then Taylor expand the symbolic functions.
%
% The original dynamics are given by:
% d/dt(x_US) = f_US(x_US) + gU_US(x_US)*u_US + gD_US(x_US)*d_US
% 
% With: x_US in xRange
%       u_US in uRange
%       d_US in dRange
%
% The scaled and approximated dynamics are given by:
% d/dt(x_SC) ~ f_SC(x_SC) + gU_SC(x_SC)*u_SC + gD_SC(x_SC)*d_SC
%
% With: x_SC in unit box
%       u_SC in unit box
%       d_SC in unit box
%
% Affine transformations are used to go from scaled to unscaled variables

if isempty(gU_US)
    gU_US = zeros(length(f_US));
    uRange = zeros(1,2);
end

if isempty(gD_US)
    gD_US = zeros(length(f_US),1);
    dRange = zeros(1,2);
end

%% Proceed through modes sequentially
%% Scale dynamics
% Get domain size and offset
xSize = diff(xRange,[],2)/2;
xMid  = mean(xRange,2);

uSize = diff(uRange,[],2)/2;
uMid = mean(uRange,2);

dSize = diff(dRange,[],2)/2;
dMid = mean(dRange,2);

% Initialize scaled x variable
x_SC_Sym = sym('xSC',[length(x_US_Sym),1]);

% Get scaling functions
x_US_Expr = diag(xSize)*x_SC_Sym + xMid; % x unscaled in terms of x scaled
x_SC_Expr = diag(1./xSize)*(x_US_Sym - xMid); % x scaled in terms of x unscaled

x_US_Expr_Poly = diag(xSize)*x_SC_Poly + xMid; % x unscaled in terms of x scaled

x_US_Fun = matlabFunction(x_US_Expr,'vars',{x_SC_Sym});
x_SC_Fun = matlabFunction(x_SC_Expr,'vars',{x_US_Sym});

u_US_Fun = @(u_SC) diag(uSize)*u_SC + uMid;
u_SC_Fun = @(u_US) diag(1./uSize)*(u_US - uMid);

d_US_Fun = @(d_SC) diag(dSize)*d_SC + dMid;
d_SC_Fun = @(d_US) diag(1./dSize)*(d_US - dMid);
    
% Get scaled f
f_SC = diag(1./xSize) * subs(f_US + gU_US*uMid + gD_US*dMid, x_US_Sym, x_US_Expr);

% Get scaled gU
gU_SC =  diag(1./xSize) * subs(gU_US*diag(uSize), x_US_Sym, x_US_Expr);

% Get scaled gD
gD_SC =  diag(1./xSize) * subs(gD_US*diag(dSize), x_US_Sym, x_US_Expr);

    
    %% Taylor approximate dynamics
% Taylor approximate f
f = msspoly(zeros(size(f_SC)));
for i = 1:length(f(:))
    f(i) = msstaylor(f_SC(i),x_SC_Sym,x_SC_Poly,0*xMid,degree+1);
end

% Taylor approximate gU
gU = msspoly(zeros(size(gU_SC)));
for i = 1:length(gU(:))
    gU(i) = msstaylor(gU_SC(i),x_SC_Sym,x_SC_Poly,0*xMid,degree+1);
end

% Taylor approximate gD
gD = msspoly(zeros(size(gD_SC)));
for i = 1:length(gD(:))
    gD(i) = msstaylor(gD_SC(i),x_SC_Sym,x_SC_Poly,0*xMid,degree+1);
end

scaleParam.x_US_Sym = x_US_Sym;
scaleParam.x_SC_Sym = x_SC_Sym;
scaleParam.x_US_Expr = x_US_Expr;
scaleParam.x_SC_Expr = x_SC_Expr;
scaleParam.x_US_Expr_Poly = x_US_Expr_Poly;

scaleParam.x_US_Fun = x_US_Fun;
scaleParam.x_SC_Fun = x_SC_Fun;
scaleParam.u_US_Fun = u_US_Fun;
scaleParam.u_SC_Fun = u_SC_Fun;
scaleParam.d_US_Fun = d_US_Fun;
scaleParam.d_SC_Fun = d_SC_Fun;

scaleParam.xSize = xSize;
scaleParam.xMid = xMid;
scaleParam.uSize = uSize;
scaleParam.uMid = uMid;
scaleParam.dSize = dSize;
scaleParam.dMid = dMid;

end