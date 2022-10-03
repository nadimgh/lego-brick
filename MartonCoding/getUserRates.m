function rate = getUserRates(sum_rate, rate_cell, scheme)
%get the rate of the users given a sum-rate and a particulat coding scheme
%  rate: rates of constituent codes at the power level that achieves (sum_rate + some rate backoff parameter) for the
%  considered coding scheme

rate_dict = (0:256)/256;
R = size(rate_cell{1}.H, 1);
rate_backoff = 1/16;

if strcmp(scheme, 'marton_opt')
    rate_marton_opt_precoding = rate_cell{1}.rate_marton_opt_precoding;
    rate_marton_opt_precoding_user = rate_cell{2};
    [~, idx] = min(abs(rate_marton_opt_precoding - sum_rate - R*rate_backoff));
    rate_unnormalized = rate_marton_opt_precoding_user{idx};
    rate = zeros(R,2);
    for i_user = 1:R
        rate(i_user,1) = rate_dict(find(rate_dict <= rate_unnormalized(i_user,1), 1, 'last')) - rate_backoff;
        rate(i_user,2) = rate_dict(find(rate_dict >= rate_unnormalized(i_user,2), 1));
    end
elseif strcmp(scheme, 'marton')
    rate_marton = rate_cell{1}.rate_marton;
    rate_marton_user = rate_cell{2};
    [~, idx] = min(abs(rate_marton - sum_rate - R*rate_backoff));
    rate_unnormalized = rate_marton_user{idx};
    rate = zeros(R,2);
    for i_user = 1:R
        rate(i_user,1) = rate_dict(find(rate_dict <= rate_unnormalized(i_user,1), 1, 'last')) - rate_backoff;
        rate(i_user,2) = rate_dict(find(rate_dict >= rate_unnormalized(i_user,2), 1));
    end
elseif strcmp(scheme, 'symmetric_opt')
    rate_symmetric_opt_precoding = rate_cell{1}.rate_symmetric_opt_precoding;
    rate_symmetric_opt_precoding_user = rate_cell{2};
    [~, idx] = min(abs(rate_symmetric_opt_precoding - sum_rate - R*rate_backoff));
    rate_unnormalized = rate_symmetric_opt_precoding_user{idx};
    rate = zeros(R,1);
    for i_user = 1:R
        rate(i_user) = rate_dict(find(rate_dict <= rate_unnormalized(i_user), 1, 'last')) - rate_backoff;
    end
elseif strcmp(scheme, 'symmetric')
    rate_symmetric = rate_cell{1}.rate_symmetric;
    rate_symmetric_user = rate_cell{2};
    [~, idx] = min(abs(rate_symmetric - sum_rate));
    rate_unnormalized = rate_symmetric_user{idx};
    rate = zeros(R,1);
    for i_user = 1:R
        rate(i_user) = rate_dict(find(rate_dict <= rate_unnormalized(i_user), 1, 'last'));
    end
elseif strcmp(scheme, 'mmse')
    rate_mmse = rate_cell{1}.rate_mmse;
    rate_mmse_user = rate_cell{2};
    [~, idx] = min(abs(rate_mmse - sum_rate));
    rate_unnormalized = rate_mmse_user{idx};
    rate = zeros(R,1);
    for i_user = 1:R
        rate(i_user) = rate_dict(find(rate_dict <= rate_unnormalized(i_user), 1, 'last'));
    end
elseif strcmp(scheme, 'zf')
    rate_zf = rate_cell{1}.rate_zf;
    rate_zf_user = rate_cell{2};
    [~, idx] = min(abs(rate_zf - sum_rate));
    rate_unnormalized = rate_zf_user{idx};
    rate = zeros(R,1);
    for i_user = 1:R
        rate(i_user) = rate_dict(find(rate_dict <= rate_unnormalized(i_user), 1, 'last'));
    end
elseif strcmp(scheme, 'time_division')
    rate_unnormalized = min(sum_rate, 1);
    rate = rate_dict(find(rate_dict <= rate_unnormalized, 1, 'last'));
else
    error('The entered coding scheme string is not valid.');
end

end

