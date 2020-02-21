function [y_SC] = unscale(y_US, y_US_Expr, y_SC_Sym)
% Unscale a variable that was scaled by scaleTaylorDynamicsHybrid

if ~iscell(y_US)
    y_US = {y_US};
    y_US_Expr = {y_US_Expr};
    y_SC_Sym = {y_SC_Sym};
end
    
y_SC = cell(size(y_US));

for i = 1:length(y_US)
    y_SC{i} = double(subs(y_US_Expr{i}', num2cell(y_SC_Sym{i})', num2cell(y_US{i},1)));
end

end