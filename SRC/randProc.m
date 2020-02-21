function [ signal ] = randProc( t, size, freq )
%RANDPROC Simulates a function in L^\infty: R -> [0,1]
rng(freq*t);
val = rand(size(1),size(2));
end

