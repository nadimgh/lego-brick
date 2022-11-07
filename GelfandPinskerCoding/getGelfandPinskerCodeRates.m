function [rate, target_dist, idx_target_snr] = getGelfandPinskerCodeRates(overall_rate, rate_backoff, p_s, g, P_dB)
%get the rates of the constituent point-to-point channel codes for the
%Gelfand-Pinkser coding scheme, where:
%  overall_rate: overall rate of the coding scheme
%  p_s: 1x2 vector for the channel state distribution p(s) (has the form [1-theta theta])
%  g: channel parameter (Y = X+gS+Z)
%  rate_backoff: backoff parameter away from the theoretical limits
%  P_dB: power constraint (in dB)

rate_dict = (0:256)/256;
rate1_ub = zeros(1, length(P_dB));
rate2_lb = zeros(1, length(P_dB));

theta = p_s(2);
matfile_capacity = sprintf('AchievableRates_GelfandPinsker_theta%.2f_g%.2f.mat', theta, g);
if exist(matfile_capacity, 'file')
    res_capacity = load(matfile_capacity);
    target_dist_all = res_capacity.target_dist;
    P_dB_all = res_capacity.P_dB;
    capacity_all = res_capacity.capacity;
    [~, idx] = ismember(P_dB, P_dB_all);
    target_dist = target_dist_all(idx,:);
    capacity = capacity_all(idx);
else
    [target_dist, target_dist_bsc, capacity, capacity_bsc, capacity_symmetric] = computeCapacityAchievingDistributionGaussian(p_s, g, P_dB);
    save(matfile_capacity, 'p_s', 'g', 'P_dB', 'target_dist', 'target_dist_bsc', 'capacity', 'capacity_bsc', 'capacity_symmetric');
end

p_x = [target_dist(:,1)+target_dist(:,3) target_dist(:,2)+target_dist(:,4)];
for i = 1:length(P_dB)
    rate2_lb(i) = computeMutualInformation(target_dist(i,:)) + 1 - h2(p_x(i,2)); % 1-H(X|S)
    rate1_ub(i) = capacity(i) + rate2_lb(i); % 1-H(X|Y)
end

[~, idx_target_snr] = min(abs(capacity - overall_rate - rate_backoff));
rate = zeros(1,2);
rate(1) = rate_dict(find(rate_dict <= rate1_ub(idx_target_snr), 1, 'last')) - rate_backoff;
rate(2) = rate_dict(find(rate_dict >= rate2_lb(idx_target_snr), 1));


end

