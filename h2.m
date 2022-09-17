function res = h2(alpha)
% binary entropy function

if ( any(alpha < 0, 'all') || any(alpha > 1, 'all') )
    error('Alpha should be in [0,1].');
end

res = zeros(size(alpha));
idx = (alpha ~= 0) & (alpha ~= 1);
res(idx) = -alpha(idx).*log2(alpha(idx)) - (1-alpha(idx)).*log2(1-alpha(idx));

end
