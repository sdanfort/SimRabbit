function pmMat = pmCombo(n)
    pmMat = zeros(n,2^n);
    for i = 1:n
        pmRow = [ones(1,2^(i-1)),-ones(1,2^(i-1))];
        pmMat(i,:) = repmat(pmRow,1,2^n/(2^(i)));
    end
end