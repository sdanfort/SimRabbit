function [ xMerge ] = mergeModes( x, m )
%MERGEMODES Summary of this function goes here
%   Detailed explanation goes here
xMerge = zeros( size( x{1} ) );
for i = 1:length(x)
    xMerge( m==i, : ) = x{i}( m==i, : );
end

end

