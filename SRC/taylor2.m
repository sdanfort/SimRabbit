function [poly,coeff,exps] = taylor2(expr, symIn, x0, order)
symIn = symIn(:).';
x0 = x0(:).';

n = length(symIn);
dMax = order-1;

% Constant term of Taylor expansion
poly = subs(expr,symIn,x0);
exps = zeros(1,n);
coeff = double(poly);

derivsOld = expr;
expsOld = exps;

for d = 1:dMax
    fprintf('Working on degree %i/%i...\n',d,dMax)
    expsNew = addDegree(expsOld);
    derivsNew = sym(zeros(size(expsOld,1),1));
    
    for i = 1:size(expsNew,1)
        fprintf('%i/%i, ',i,size(expsNew,1))
        exp_ = expsNew(i,:);
        % Get term from previous degree
        [maxExp,indDiff] = max(exp_);
        expFind_ = exp_;
        expFind_(indDiff) = maxExp-1;
        
        indOld = find(prod((expsOld == expFind_).'),1,'first');
        
        derivsNew(i) = diff(derivsOld(indOld),symIn(indDiff));
        coeff_ = (subs(derivsNew(i),symIn,x0))/prod(factorial(exp_));
        
        coeff = [coeff;double(coeff_)];
        poly = poly + coeff_*prod((symIn - x0).^exp_);
    end
    fprintf('\n')
    
    exps = [exps;expsNew];
    
    expsOld = expsNew;
    derivsOld = derivsNew;
end


end

function expsNew = addDegree(expsOld)
    l = size(expsOld,1);
    n = size(expsOld,2);
    
    indVec = repmat(1:l,n,1);
    expsNew = expsOld(indVec,:) + repmat(eye(n),l,1);
    
%     expsNew = zeros(n*size(expsOld,1),n);
%     % Add 1 to each element of each row of exps Old
%     for i = 1:l
%         expsNew((i-1)*n + (1:n),:) = repmat(expsOld(i,:),n,1);
%         for j = 1:n
%             expsNew((i-1)*n + j,j) = expsNew((i-1)*n + j,j)+1;
%         end
%     end

    % Remove duplicates
    expsNew = unique(expsNew,'rows');
end