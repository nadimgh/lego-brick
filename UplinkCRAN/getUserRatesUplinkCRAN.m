function rate = getUserRatesUplinkCRAN(sum_rate, rate_backoff, constraint_backoff, rate_cell)
%get the rate of the users given a sum-rate
%  rate: rates of constituent codes at the power level that achieves sum_rate
%  rate_backoff: backoff parameter away from the theoretical limits of error correction codes at the user side
%  constraint_backoff: backoff parameter away from the theoretical limits error correction codes at the relays

rate_dict = (0:256)/256;
[T, R] = size(rate_cell{1}.H); % R-user, T-relay

max_sum_rate = rate_cell{1}.max_sum_rate;
rate_user = rate_cell{2};
[~, idx] = min(abs(max_sum_rate - sum_rate - R*rate_backoff));
rate_unnormalized = rate_user{idx};

rate = zeros(R+T,2);
for i_user = 1:R
    rate(i_user,1) = rate_dict(find(rate_dict <= rate_unnormalized(i_user,1), 1, 'last')) - rate_backoff;
    rate(i_user,2) = rate_dict(find(rate_dict >= rate_unnormalized(i_user,2), 1));
end
for i_relay = R+1:R+T
    temp = rate_dict(find(rate_dict >= rate_unnormalized(i_relay,2), 1));
    rate(i_relay,2) = max(temp - constraint_backoff, 0);
    rate(i_relay,1) = rate_dict(find(rate_dict <= rate_unnormalized(i_relay,1), 1, 'last'));
end

end

