function [rate_cell, dist_struct] = getMaximumAchievableUserRatesUplinkCRAN(H, P_dB, C_vec)
%get the maximum achievable rates and capacity achieving distributions for the uplink C-RAN coding
%scheme, where:
%  rate_cell: 1x2 cell containing parameters of the maximum achievable rates 
%  rate_cell{1}: struct containing H, P_dB, maximum achievable sum-rate
%  rate_cell{2}: 1xlength(P_dB) cell each containing the maximum achievable user rates
%  dist_struct: struct containing the capacity-achieving distributions used by each
%  constituent code

g = H(1,2);
[T,R] = size(H); % R-user, T-relay
P = 10.^(P_dB/10);

matfile = sprintf('AchievableRates_UplinkCRAN_g%.2f_C%.2f.mat', g, C_vec(1));

if exist(matfile, 'file')
    res = load(matfile);
else
    error('Matfile does not exist.');
end

P_dB_all = res.P_dB;
[~, idx] = ismember(P_dB, P_dB_all);

input_dist_all = res.input_dist;
joint_input_quantization_dist_all = res.joint_input_quantization_dist;
joint_output_quantization_dist_all = res.joint_output_quantization_dist;
sum_rate_all = res.sum_rate;

input_dist = input_dist_all(idx,:);
joint_input_quantization_dist = joint_input_quantization_dist_all(idx,:);
joint_output_quantization_dist = joint_output_quantization_dist_all(idx,:);
sum_rate = sum_rate_all(idx);

rate_cell = cell(1,2);
rate_user = cell(1,length(P));
p_x1 = zeros(length(P), 2);
p_x2 = zeros(length(P), 2);
p_u1u2x1 = zeros(length(P), 8); % attention to the order
p_x1u1u2x2 = zeros(length(P), 16);
p_y1u1 = zeros(length(P), 4);
p_y2u2 = zeros(length(P), 4);
p_u1 = zeros(length(P), 2);
p_u2 = zeros(length(P), 2);
p_u1u2 = zeros(length(P), 4);
for i_p = 1:length(P)
    % names of variables are indicative of the distributions that are computed ((U1,U2) --> (Y1hat,Y2hat))
    % E.g. p_x1u1u2x2 = [P(X1=0,Y1hat=0,Y2hat=0,X2=0) P(X1=0,Y1hat=0,Y2hat=0,X2=1) P(X1=0,Y1hat=0,Y2hat=1,X2=0) P(X1=0,Y1hat=0,Y2hat=1,X2=1) ...
    %                    P(X1=0,Y1hat=1,Y2hat=0,X2=0) P(X1=0,Y1hat=1,Y2hat=0,X2=1) P(X1=0,Y1hat=1,Y2hat=1,X2=0) P(X1=0,Y1hat=1,Y2hat=1,X2=1) ...
    %                    P(X1=1,Y1hat=0,Y2hat=0,X2=0) P(X1=1,Y1hat=0,Y2hat=0,X2=1) P(X1=1,Y1hat=0,Y2hat=1,X2=0) P(X1=1,Y1hat=0,Y2hat=1,X2=1) ...
    %                    P(X1=1,Y1hat=1,Y2hat=0,X2=0) P(X1=1,Y1hat=1,Y2hat=0,X2=1) P(X1=1,Y1hat=1,Y2hat=1,X2=0) P(X1=1,Y1hat=1,Y2hat=1,X2=1) ]
    
    p_x1(i_p,:) = sum(reshape(input_dist(i_p,:), 2, []));
    
    p_x2(i_p,:) = transpose(sum(reshape(input_dist(i_p,:), 2, []), 2));
    
    p_u1u2x1(i_p,:) = [sum(joint_input_quantization_dist(i_p,[1 5])) sum(joint_input_quantization_dist(i_p,[9 13])) ...
        sum(joint_input_quantization_dist(i_p,[2 6])) sum(joint_input_quantization_dist(i_p,[10 14])) ...
        sum(joint_input_quantization_dist(i_p,[3 7])) sum(joint_input_quantization_dist(i_p,[11 15])) ...
        sum(joint_input_quantization_dist(i_p,[4 8])) sum(joint_input_quantization_dist(i_p,[12 16]))];
    
    p_x1u1u2x2(i_p,:) = joint_input_quantization_dist(i_p,[1 5 2 6 3 7 4 8 9 13 10 14 11 15 12 16]);
    
    p_y1u1(i_p,:) = [sum(joint_output_quantization_dist(i_p,[1 2 5 6])) sum(joint_output_quantization_dist(i_p,[3 4 7 8])) ...
        sum(joint_output_quantization_dist(i_p,[9 10 13 14])) sum(joint_output_quantization_dist(i_p,[11 12 15 16]))];
    
    p_y2u2(i_p,:) = [sum(joint_output_quantization_dist(i_p,[1 3 9 11])) sum(joint_output_quantization_dist(i_p,[2 4 10 12])) ...
        sum(joint_output_quantization_dist(i_p,[5 7 13 15])) sum(joint_output_quantization_dist(i_p,[6 8 14 16]))];
    
    p_u1(i_p,:) = [sum(joint_output_quantization_dist(i_p,[1 2 5 6 9 10 13 14])) sum(joint_output_quantization_dist(i_p,[3 4 7 8 11 12 15 16]))];
    
    p_u2(i_p,:) = [sum(joint_output_quantization_dist(i_p,[1 3 5 7 9 11 13 15])) sum(joint_output_quantization_dist(i_p,[2 4 6 8 10 12 14 16]))];
    
    p_u1u2(i_p,:) = [sum(joint_output_quantization_dist(i_p,[1 5 9 13])) sum(joint_output_quantization_dist(i_p,[2 6 10 14])) ...
        sum(joint_output_quantization_dist(i_p,[3 7 11 15])) sum(joint_output_quantization_dist(i_p,[4 8 12 16]))];
    
    rate_user{i_p} = zeros(R+T,2);
    rate_user{i_p}(1,2) = 1 - h2(p_x1(i_p,2)); % 1 - H(X1)
    rate_user{i_p}(2,2) = 1 - h2(p_x2(i_p,2)); % 1 - H(X2)
    rate_user{i_p}(1,1) = rate_user{i_p}(1,2) + computeMutualInformation(p_u1u2x1(i_p,:)); % 1 - H(X1 | Y1hat,Y2hat)
    rate_user{i_p}(2,1) = rate_user{i_p}(2,2) + computeMutualInformation(p_x1u1u2x2(i_p,:)); % 1 - H(X2 | X1,Y1hat,Y2hat)
    rate_user{i_p}(3,2) = 1 - h2(p_u1(i_p,2)); % 1 - H(Y1hat)
    rate_user{i_p}(4,2) = computeMutualInformation(p_u1u2(i_p,:)) + 1 - h2(p_u2(i_p,2)); % 1 - H(Y2hat | Y1hat)
    rate_user{i_p}(3,1) = rate_user{i_p}(3,2) + computeMutualInformation(p_y1u1(i_p,:)); % 1 - H(Y1hat | Y1)
    rate_user{i_p}(4,1) = 1 - h2(p_u2(i_p,2)) + computeMutualInformation(p_y2u2(i_p,:)); % 1 - H(Y2hat | Y2)
end

rate_cell{1} = struct('H', H, 'P_dB', P_dB, 'max_sum_rate', sum_rate);
rate_cell{2} = rate_user;
dist_struct = struct('asym_encoder1', p_x1, 'asym_encoder2', p_x2, 'asym_decoder1', p_u1u2x1, 'asym_decoder2', p_x1u1u2x2, ...
    'lossy_encoder', p_y1u1, 'lossy_decoder', p_u1, 'wz_encoder', p_y2u2, 'wz_decoder', p_u1u2);

end