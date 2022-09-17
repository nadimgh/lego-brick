function P_opt = getOptimalZeroForcingPowerAllocation(H, P)
%get the optimal power allocation for zero forcing under an average power
%constraint P using a water filling algorithm. P_opt is a vector of length R, 
%where R is the number of users in the broadcast channel with channel matrix H.

[R, ~] = size(H);

gamma = diag(inv(H*H'));
[~, idx] = sort(gamma);

i = 1;
water_filled = 0;
mu = gamma(idx(1));
while(i<R && i*(gamma(idx(i+1)) - gamma(idx(i))) < P-water_filled)
    mu = mu + gamma(idx(i+1)) - gamma(idx(i));
    water_filled = water_filled + i*(gamma(idx(i+1)) - gamma(idx(i)));
    i = i+1;
end

mu = mu + (P-water_filled)/i;
P_opt = mu./gamma - 1;
P_opt(P_opt < 0) = 0;


end

