function [f, gU, gD, S, R, scaleParam] = scaleTaylorDynamicsHybrid(x_US_Sym, x_SC_Poly, f_US, gU_US, gD_US, S_US, R_US, xRange, uRange, dRange, degree, x0_US)
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

%% Initialize variables
nModes = length(x_US_Sym);

f = cell(nModes,1);
gU = cell(nModes,1);
gD = cell(nModes,1);
S = cell(nModes,nModes);
R = cell(nModes,nModes);

x_SC_Sym = cell(nModes,1);
x_US_Expr = cell(nModes,1);
x_SC_Expr = cell(nModes,1);
x_US_Expr_Poly = cell(nModes,1);
x_US_Fun = cell(nModes,1);
x_SC_Fun = cell(nModes,1);
u_US_Fun = cell(nModes,1);
u_SC_Fun = cell(nModes,1);
d_US_Fun = cell(nModes,1);
d_SC_Fun = cell(nModes,1);
xSize = cell(nModes,1);
xMid = cell(nModes,1);
x0_SC = cell(nModes,1);
uSize = cell(nModes,1);
uMid = cell(nModes,1);
dSize = cell(nModes,1);
dMid = cell(nModes,1);

if ~isstruct(degree)
    d_ = degree;
    degree = [];
    degree.f = d_;
    degree.gU = d_;
    degree.gD = d_;
    degree.S = d_;
    degree.R = d_;
end
    

%% Proceed through modes sequentially
for m = 1:nModes
    fprintf('Scaling dyamics, working on mode %i/%i\n',m,nModes)
    
    %% Scale dynamics
    % Get domain size and offset
    xSize{m} = diff(xRange{m},[],2)/2;
    xMid{m}  = mean(xRange{m},2);

    uSize{m} = diff(uRange{m},[],2)/2;
    uMid{m} = mean(uRange{m},2);

    dSize{m} = diff(dRange{m},[],2)/2;
    dMid{m} = mean(dRange{m},2);

    % Initialize scaled x variable
    x_SC_Sym{m} = sym('xSC',[length(x_US_Sym{m}),1]);
    
    % Get scaling functions
    x_US_Expr{m} = diag(xSize{m})*x_SC_Sym{m} + xMid{m}; % x unscaled in terms of x scaled
    x_SC_Expr{m} = diag(1./xSize{m})*(x_US_Sym{m} - xMid{m}); % x scaled in terms of x unscaled
    
    x_US_Expr_Poly{m} = diag(xSize{m})*x_SC_Poly{m} + xMid{m}; % x unscaled in terms of x scaled

    x_US_Fun{m} = matlabFunction(x_US_Expr{m},'vars',x_SC_Sym(m));
    x_SC_Fun{m} = matlabFunction(x_SC_Expr{m},'vars',x_US_Sym(m));
    
    u_US_Fun{m} = @(u_SC) diag(uSize{m})*u_SC + uMid{m};
    u_SC_Fun{m} = @(u_US) diag(1./uSize{m})*(u_US - uMid{m});
    
    d_US_Fun{m} = @(d_SC) diag(dSize{m})*d_SC + dMid{m};
    d_SC_Fun{m} = @(d_US) diag(1./dSize{m})*(d_US - dMid{m});
    
    % Get scaled f, gD, GU
    f_SC = diag(1./xSize{m}) * subs(f_US{m}, x_US_Sym{m}, x_US_Expr{m});
    gU_SC = [];
    gD_SC = [];
    
    if ~isempty(gU_US{m})
        % Get scaled gU
        gU_SC =  diag(1./xSize{m}) * subs(gU_US{m}*diag(uSize{m}), x_US_Sym{m}, x_US_Expr{m});
        f_SC = f_SC + diag(1./xSize{m}) * subs(gU_US{m}*uMid{m}, x_US_Sym{m}, x_US_Expr{m});
    end
    
    if ~isempty(gD_US{m})
        % Get scaled gD
        gD_SC =  diag(1./xSize{m}) * subs(gD_US{m}*diag(dSize{m}), x_US_Sym{m}, x_US_Expr{m});
        f_SC = f_SC + diag(1./xSize{m}) * subs(gD_US{m}*dMid{m}, x_US_Sym{m}, x_US_Expr{m});
    end
    
    %% Taylor approximate dynamics
    if nargin >= 12
        x0_SC{m} = x_SC_Fun{m}(x0_US{m});
    else
        x0_SC{m} = 0*xMid{m};
    end
        
    % Taylor approximate f
    f{m} = msspoly(zeros(size(f_SC)));
    for i = 1:length(f{m}(:))
        f{m}(i) = msstaylor(f_SC(i),x_SC_Sym{m},x_SC_Poly{m},x0_SC{m},degree.f+1);
    end

    % Taylor approximate gU
    gU{m} = msspoly(zeros(size(gU_SC)));
    for i = 1:length(gU{m}(:))
        gU{m}(i) = msstaylor(gU_SC(i),x_SC_Sym{m},x_SC_Poly{m},x0_SC{m},degree.gU+1);
    end

    % Taylor approximate gD
    gD{m} = msspoly(zeros(size(gD_SC)));
    for i = 1:length(gD{m}(:))
        gD{m}(i) = msstaylor(gD_SC(i),x_SC_Sym{m},x_SC_Poly{m},x0_SC{m},degree.gD+1);
    end
end
   
%% Scale and taylor guards and reset maps
for m = 1:nModes
    fprintf('Scaling guards and resets, working on mode %i/%i\n',m,nModes)
    for mp = 1:nModes
        % Scale input of S
        S_SC = subs(S_US{m,mp},x_US_Sym{m},x_US_Expr{m});
        % Taylor S
        S{m,mp} = msspoly(zeros(size(S_SC)));
        for i = 1:length(S{m,mp}(:))
            S{m,mp}(i) = msstaylor(S_SC(i),x_SC_Sym{m},x_SC_Poly{m},x0_SC{m},degree.S+1);
        end
        
        % Scale input of R
        R_temp = subs(R_US{m,mp},x_US_Sym{m}, x_US_Expr{m});
        % Scale output of R
        R_SC = subs(x_SC_Expr{mp},x_US_Sym{mp},R_temp);
        % Taylor R
        R{m,mp} = msspoly(zeros(size(R_SC)));
        for i = 1:length(R{m,mp}(:))
            R{m,mp}(i) = msstaylor(R_SC(i),x_SC_Sym{m},x_SC_Poly{m},x0_SC{m},degree.R+1);
        end
    end
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