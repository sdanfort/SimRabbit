function [ msspolyOut ] = msstaylor( symFunIn, symIn, mssOut, x0, order )
%MSSTAYLOR Summary of this function goes here
%   Detailed explanation goes here
% 
% [~,coeff,deg] = taylor2(symFunIn, symIn, x0, order);
% 
% msspolyOut = recomp(mssOut(:),deg,coeff(:).');
% msspolyOut = subs(msspolyOut,mssOut(:),mssOut(:)-x0);


f_poly = taylor(symFunIn, symIn, x0, 'order', order);
[coeff,monom] = coeffs(f_poly,symIn);

f_poly_h = matlabFunction(double(coeff)*monom','vars',{symIn});       

msspolyOut = f_poly_h(mssOut);

end

