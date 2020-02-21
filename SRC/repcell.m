function cellArrayOut = repcell(cellArrayIn,nRow,nCol)

if ~iscell(cellArrayIn)
    cellArrayIn = {cellArrayIn};
end

h = size(cellArrayIn,1);
w = size(cellArrayIn,2);

cellArrayOut = cell(h*nRow,w*nCol);

for i = 0:nRow-1
    for j = 0:nCol-1
        cellArrayOut((i*h+1):((i+1)*h),(j*w+1):((j+1)*w)) = cellArrayIn;
    end
end
