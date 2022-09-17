function rate = getLossySourceCodeRates(target_dist, rate_backoff)
%get the code rates of a given target joint distribution and a rate backoff
%parameter

p_xhat_target = sum(reshape(target_dist, 2, []), 2);
rate_dict = (0:256)/256;

rate = zeros(1,2);
rate2_bound = 1 - h2(p_xhat_target(2));
rate1_bound = computeMutualInformation(target_dist) + 1 - h2(p_xhat_target(2));
rate(2) = max(0, rate_dict(find(rate_dict <= rate2_bound, 1, 'last')) - rate_backoff);
rate(1) = rate_dict(find(rate_dict >= rate1_bound, 1));

end

