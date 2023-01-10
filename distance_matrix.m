
function [ matrix ] = distance_matrix(frag_1, tau)
% returns the distance matrix between a spike train fragment
% using the van Rossum metric
% frag1  are all the fragments of spike train 
% tau is a timescale expressing the precision of spike times

import metrics.*

n = length(frag_1);       % two fragments should have the same length
matrix = zeros(n, n);    % initialize matrix

for i = 1:n
    for j = 1:n
        interval_1 = cell2mat(frag_1(i));
        interval_2 = cell2mat(frag_1(j));
        distance = square_term(interval_1, tau) + square_term(interval_2, tau) - 2*cross_term(interval_1, interval_2, tau);
        matrix(i,j) = distance;
        matrix(j,i) = distance;

    end
end
end