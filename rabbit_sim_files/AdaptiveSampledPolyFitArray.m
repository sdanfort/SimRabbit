function coeffArray = AdaptiveSampledPolyFitArray(int_,monomsInit_,polyInit_,monomsTest_,polyTest_,boundSign,eps_,BLP,cLP)
if nargin < 7
    eps_ = 0.05;
end
if nargin < 9
    BLP = [];
    cLP = [];
end

f = boundSign*int_;
A = -boundSign*monomsInit_.';
b = -boundSign*polyInit_-eps_;

options = optimoptions('linprog','display','off');

bounded = 0;
while ~bounded
    [coeffs_,~,exitFlag] = linprog(f,A,b,BLP,cLP,[],[],options);
    
    if exitFlag == -3
        warning('Problem unbounded')
    end
    
    PolyTestUB_ = monomsTest_.'*coeffs_;
    
%     indAdd = (boundSign*(PolyTestUB_-polyTest_)<eps_/2);
%     fprintf('Number of violations: %i,  ', sum(indAdd))
%     fprintf('Max violation: %4.6f,  ', max(boundSign*(polyTest_-PolyTestUB_)))
%     fprintf('New problem size: %i\n', sum(indAdd) + length(b))

indAdd = 0;
    
    if sum(indAdd) == 0
        bounded = 1;
    else
        A = [A;
             -boundSign*monomsTest_(:,indAdd).'];
        b = [b;
             -boundSign*polyTest_(indAdd) - eps_];
    end
end

coeffArray = coeffs_;

end