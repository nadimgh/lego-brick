function rate_cell = getMaximumAchievableUserRates(H, P_dB, scheme)
%get the maximum achievable rates by each user for one of the coding
%schemes that we consider, where:
%  rate_cell: 1x2 cell (or 1x3 cell if precoding matrices are to be returned) containing parameters of the maximum achievable rates 
%  (rate_cell{1} is a struct containing H, P_dB, maximum achievable sum-rate of the coding scheme and the capacity-achieving distribution, if needed, 
%  rate_cell{2} is a 1xlength(P_dB) cell each containing the maximum achievable user rates, 
%  rate_cell{3} is a 1xlength(P_dB) cell each containing the optimal precoding matrix to be used at each power level)

g = H(1,2);
[R, T] = size(H);
P = 10.^(P_dB/10);
matfile = sprintf('AchievableRates_BC_g%.2f.mat', g);

if exist(matfile, 'file')
    res = load(matfile);
else
    error('Matfile does not exist.');
end

P_dB_all = res.P_dB;
[~, idx] = ismember(P_dB, P_dB_all);

if strcmp(scheme, 'marton_opt')
    rate_cell = cell(1,3);
    input_dist_marton_opt_precoding_all = res.input_dist_marton_opt_precoding;
    W_opt_all = res.W_opt;
    input_dist_marton_opt_precoding = zeros(length(P), 2^R);
    W_opt = cell(1, length(P));
    rate_marton_opt_precoding_user = cell(1, length(P));
    rate_marton_opt_precoding = zeros(1, length(P));
    for i_p = 1:length(P)
        input_dist_marton_opt_precoding(i_p,:) = input_dist_marton_opt_precoding_all(idx(i_p),:);
        W_opt{i_p} = W_opt_all{idx(i_p)};
        
        rate_marton_opt_precoding_user{i_p} = zeros(R,2);
        p_x1 = sum(reshape(input_dist_marton_opt_precoding(i_p,:), 2, []));
        p_x2 = transpose(sum(reshape(input_dist_marton_opt_precoding(i_p,:), 2, []), 2));
        [rate_marton_opt_precoding(i_p), rate_user] = computeCapacityMartonSpecificInputDist(H, P(i_p), W_opt{i_p}, input_dist_marton_opt_precoding(i_p,:));
        rate_marton_opt_precoding_user{i_p}(1,2) = 1 - h2(p_x1(2)); % 1-H(X1)
        rate_marton_opt_precoding_user{i_p}(2,2) = computeMutualInformation(input_dist_marton_opt_precoding(i_p,:)) + 1 - h2(p_x2(2)); % 1-H(X2|X1)
        rate_marton_opt_precoding_user{i_p}(1,1) = rate_marton_opt_precoding_user{i_p}(1,2) + abs(rate_user(1)); % 1-H(X1|Y1)
        rate_marton_opt_precoding_user{i_p}(2,1) = rate_marton_opt_precoding_user{i_p}(2,2) + abs(rate_user(2)); % 1-H(X2|Y2)
    end
    rate_cell{1} = struct('H', H, 'P_dB', P_dB, 'rate_marton_opt_precoding', rate_marton_opt_precoding, 'input_dist_marton_opt_precoding', input_dist_marton_opt_precoding);
    rate_cell{2} = rate_marton_opt_precoding_user;
    rate_cell{3} = W_opt;
elseif strcmp(scheme, 'marton')
    rate_cell = cell(1,2);
    input_dist_marton_all = res.input_dist_marton;
    input_dist_marton = zeros(length(P), 2^R);
    rate_marton_user = cell(1, length(P));
    rate_marton = zeros(1, length(P));
    for i_p = 1:length(P)
        input_dist_marton(i_p,:) = input_dist_marton_all(idx(i_p),:);
        rate_marton_user{i_p} = zeros(R,2);
        p_x1 = sum(reshape(input_dist_marton(i_p,:), 2, []));
        p_x2 = transpose(sum(reshape(input_dist_marton(i_p,:), 2, []), 2));
        [rate_marton(i_p), rate_user] = computeCapacityMartonSpecificInputDist(H, P(i_p), eye(R), input_dist_marton(i_p,:));
        rate_marton_user{i_p}(1,2) = 1 - h2(p_x1(2));% 1-H(X1)
        rate_marton_user{i_p}(2,2) = computeMutualInformation(input_dist_marton(i_p,:)) + 1 - h2(p_x2(2)); % 1-H(X2|X1)
        rate_marton_user{i_p}(1,1) = rate_marton_user{i_p}(1,2) + abs(rate_user(1)); % 1-H(X1|Y1)
        rate_marton_user{i_p}(2,1) = rate_marton_user{i_p}(2,2) + abs(rate_user(2)); % 1-H(X2|Y2)
    end
    rate_cell{1} = struct('H', H, 'P_dB', P_dB, 'rate_marton', rate_marton, 'input_dist_marton', input_dist_marton);
    rate_cell{2} = rate_marton_user;
elseif strcmp(scheme, 'symmetric_opt')
    rate_cell = cell(1,3);
    W_opt_symmetric_all = res.W_opt_symmetric;
    W_opt_symmetric = cell(1, length(P));
    rate_symmetric_opt_precoding_user = cell(1, length(P));
    rate_symmetric_opt_precoding = zeros(1, length(P));
    for i_p = 1:length(P)
        W_opt_symmetric{i_p} = W_opt_symmetric_all{idx(i_p)};
        
        [rate_symmetric_opt_precoding(i_p), rate_user] = computeCapacityPrecoding(H, P(i_p), W_opt_symmetric{i_p});
        rate_symmetric_opt_precoding_user{i_p} = abs(rate_user');
    end
    rate_cell{1} = struct('H', H, 'P_dB', P_dB, 'rate_symmetric_opt_precoding', rate_symmetric_opt_precoding);
    rate_cell{2} = rate_symmetric_opt_precoding_user;
    rate_cell{3} = W_opt_symmetric;
elseif strcmp(scheme, 'symmetric')
    rate_cell = cell(1,2);
    rate_symmetric_user = cell(1, length(P));
    rate_symmetric = zeros(1, length(P));
    for i_p = 1:length(P)
        [rate_symmetric(i_p), rate_user] = computeCapacityPrecoding(H, P(i_p), eye(R));
        rate_symmetric_user{i_p} = abs(rate_user');
    end
    rate_cell{1} = struct('H', H, 'P_dB', P_dB, 'rate_symmetric', rate_symmetric);
    rate_cell{2} = rate_symmetric_user;
elseif strcmp(scheme, 'mmse')
    rate_cell = cell(1,2);
    rate_mmse_user = cell(1, length(P));
    rate_mmse = zeros(1, length(P));
    for i_p = 1:length(P)
        W_mmse = (H'*H + (R/P(i_p))*eye(T))\(H');
        [rate_mmse(i_p), rate_user] = computeCapacityPrecoding(H, P(i_p), W_mmse);
        rate_mmse_user{i_p} = abs(rate_user');
    end
    rate_cell{1} = struct('H', H, 'P_dB', P_dB, 'rate_mmse', rate_mmse);
    rate_cell{2} = rate_mmse_user;
elseif strcmp(scheme, 'zf')
    rate_cell = cell(1,2);
    rate_zf_user = cell(1, length(P));
    rate_zf = zeros(1, length(P));
    for i_p = 1:length(P)
        [rate_zf(i_p), rate_user] = computeCapacityZeroForcing(H, P(i_p));
        rate_zf_user{i_p} = abs(rate_user');
    end
    rate_cell{1} = struct('H', H, 'P_dB', P_dB, 'rate_zf', rate_zf);
    rate_cell{2} = rate_zf_user;
elseif strcmp(scheme, 'time_division')
    rate_cell = cell(1,1);
    rate_time_division = zeros(1, length(P));
    for i_p = 1:length(P)
        rate_time_division(i_p) = computeCapacityTimeDivision(H, P(i_p));
    end
    rate_cell{1} = struct('H', H, 'P_dB', P_dB, 'rate_time_division', rate_time_division);
else
    error('The entered coding scheme string is not valid.');
end


end

