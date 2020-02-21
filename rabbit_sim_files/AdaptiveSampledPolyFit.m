function polyUB = AdaptiveSampledPolyFit(monoms,int_,monomsInit_,polyInit_,monomsTest_,polyTest_,boundSign,eps_)
if nargin < 8
    eps_ = 0.05;
end
nmonom = length(monoms);

f = boundSign*int_;
A = -boundSign*monomsInit_.';
b = -boundSign*polyInit_-eps_;

options = optimoptions('linprog','display','off');

bounded = 0;
while ~bounded
    coeffs_ = linprog(f,A,b,[],[],[],[],options);
    
    PolyTestUB_ = monomsTest_.'*coeffs_(1:nmonom);
    
    indAdd = (boundSign*(PolyTestUB_-polyTest_)<eps_/2);
    fprintf('Number of violations: %i,  ', sum(indAdd))
    fprintf('Max violation: %4.6f,  ', max(boundSign*(polyTest_-PolyTestUB_)))
    fprintf('New problem size: %i\n', sum(indAdd) + length(b))
    
    if sum(indAdd) == 0
        bounded = 1;
        
        polyUB = coeffs_(1:nmonom).'*monoms;
    else
        A = [A;
             -boundSign*monomsTest_(:,indAdd).'];
        b = [b;
             -boundSign*polyTest_(indAdd) - eps_];
    end
end

end