function rate = getUserRatesDownlinkCRAN(sum_rate, rate_backoff, constraint_backoff, rate_cell, scheme)
%get the rate of the users given a sum-rate
%  rate: rates of constituent codes at the power level that achieves sum_rate
%  rate_backoff: backoff parameter away from the theoretical limits for channel codes
%  constraint_backoff: backoff parameter away from the theoretical limits for compression codes

rate_dict = (0:256)/256;
[R, T] = size(rate_cell{1}.H);

if strcmp(scheme, 'cran')
    rate_cran = rate_cell{1}.rate_cran;
    rate_cran_user = rate_cell{2};
    [~, idx] = min(abs(rate_cran - sum_rate - R*rate_backoff));
    rate_unnormalized = rate_cran_user{idx};
    rate = zeros(R+T,2);
    for i_user = 1:R
        rate(i_user,1) = rate_dict(find(rate_dict <= rate_unnormalized(i_user,1), 1, 'last')) - rate_backoff;
        rate(i_user,2) = rate_dict(find(rate_dict >= rate_unnormalized(i_user,2), 1));
    end
    for i_relay = R+1:R+T
        temp = rate_dict(find(rate_dict >= rate_unnormalized(i_relay,2), 1));
%         rate_diff = min(temp, constraint_backoff);
        rate(i_relay,2) = max(temp - constraint_backoff, 0);
        rate(i_relay,1) = rate_dict(find(rate_dict >= rate_unnormalized(i_relay,1), 1));
    end
else
    error('The entered coding scheme string is not valid.');
end

end

