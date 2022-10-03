function rate_cell = getMaximumAchievableUserRatesCRAN(H, P_dB, scheme, C_vec)
%get the maximum achievable rates by each user for one of the coding
%schemes that we consider, where:
%  rate_cell: 1x2 cell (or 1x3 cell if precoding matrices are to be returned) containing parameters of the maximum achievable rates 
%  (rate_cell{1} is a struct containing H, P_dB, maximum achievable sum-rate of the coding scheme and the capacity-achieving distribution, if needed, 
%  rate_cell{2} is a 1xlength(P_dB) cell each containing the maximum achievable user rates, 
%  rate_cell{3} is a 1xlength(P_dB) cell each containing the optimal precoding matrix to be used at each power level)

g = H(1,2);
[R, T] = size(H);
P = 10.^(P_dB/10);
matfile = sprintf('AchievableRates_DownlinkCRAN_GA_g%.2f_C%.2f.mat', g, C_vec(1));

if exist(matfile, 'file')
    res = load(matfile);
else
    error('Matfile does not exist.');
end

P_dB_all = res.P_dB;
[~, idx] = ismember(P_dB, P_dB_all);
x_all = transpose(get_tuples(R+T));
user_rate = zeros(R, length(P));

if strcmp(scheme, 'cran')
    rate_cell = cell(1,3);
    input_dist_all = res.input_dist;
    rate_cran_all = res.sum_rate;
    W_opt_all = res.W_opt;
    input_dist = zeros(length(P), 2^(R+T));
    W_opt = cell(1, length(P));
    rate_cran_user = cell(1, length(P));
    rate_cran = zeros(1, length(P));
    for i_p = 1:length(P)
        input_dist(i_p,:) = input_dist_all(idx(i_p),:);
        W_opt{i_p} = W_opt_all{idx(i_p)};
        rate_cran(i_p) = rate_cran_all(idx(i_p));
        [~, user_rate(:,i_p)] = computeCapacityCranSpecificInputDist(H, P(i_p), W_opt{i_p}, input_dist(i_p,:));
        
        idx1 = (x_all(1,:) == 1);
        p_u1 = [1-sum(input_dist(i_p,idx1)) sum(input_dist(i_p,idx1))];
        
        idx2 = (x_all(2,:) == 1);
        p_u2 = [1-sum(input_dist(i_p,idx2)) sum(input_dist(i_p,idx2))];
        
        idx3 = (x_all(3,:) == 1);
        p_x1 = [1-sum(input_dist(i_p,idx3)) sum(input_dist(i_p,idx3))];
        
        idx4 = (x_all(4,:) == 1);
        p_x2 = [1-sum(input_dist(i_p,idx4)) sum(input_dist(i_p,idx4))];
        
        p_u1u2 = sum(reshape(input_dist(i_p,:), 2^T, []));
        p_u1u2x1 = sum(reshape(input_dist(i_p,:), 2, []));
        
        rate_cran_user{i_p} = zeros(R+T,2);
        rate_cran_user{i_p}(1,2) = 1 - h2(p_u1(2));
        rate_cran_user{i_p}(2,2) = computeMutualInformation(p_u1u2) + 1 - h2(p_u2(2));
        rate_cran_user{i_p}(1,1) = rate_cran_user{i_p}(1,2) + user_rate(1,i_p);
        rate_cran_user{i_p}(2,1) = rate_cran_user{i_p}(2,2) + user_rate(2,i_p);
        rate_cran_user{i_p}(3,1) = computeMutualInformation(p_u1u2x1) + 1 - h2(p_x1(2));
        rate_cran_user{i_p}(3,2) = 1 - h2(p_x1(2));
        rate_cran_user{i_p}(4,1) = computeMutualInformation(input_dist(i_p,:)) + 1 - h2(p_x2(2));
        rate_cran_user{i_p}(4,2) = 1 - h2(p_x2(2));
    end
    rate_cell{1} = struct('H', H, 'P_dB', P_dB, 'rate_cran', rate_cran, 'input_dist', input_dist);
    rate_cell{2} = rate_cran_user;
    rate_cell{3} = W_opt;
else
    error('The entered coding scheme string is not valid.');
end


end

