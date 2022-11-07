function entropy = computeEntropy(input_dist)
%compute the entropy of the input distribution input_dist

if ( any(input_dist < 0) || any(input_dist > 1) || abs(sum(input_dist) - 1) > 1e-5 )
    error('Input_dist should be a probability distribution.');
end
input_dist = input_dist(:);

input_dist = input_dist(input_dist ~= 0);
entropy = sum( -1 * input_dist .* log2(input_dist) );

end
