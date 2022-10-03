function [input_dist, W_opt, sum_rate, W_opt_symmetric, sum_rate_symmetric] = computeCapacityMartonWithOptimalPrecoding(H, P, W_initial)
%find the optimal precoding matrix that maximizes the sum-rate
%I(X1;Y1)+I(X2;Y2)-I(X1;X2)
%for the Gaussian broadcast channel with channel gain matrix H, power
%constraint P, and input distribution input_dist

tic
if nargin < 3
    W_initial = randn(2,2);
end

[R, ~] = size(H);
alpha = 0.02:0.02:0.5;
beta = 0.02:0.02:0.5;
gamma = 0.02:0.02:0.5;
W_initial_tr = transpose(W_initial);

s_all = 1-2*transpose(get_tuples(R));

low_val = min(sqrt(P) * H * s_all, [], 'all') - 5;
upper_val = max(sqrt(P) * H * s_all, [], 'all') + 5;

input_dist_alpha = cell(length(alpha), length(beta), length(gamma));
rate_marton_alpha = zeros(length(alpha), length(beta), length(gamma));
W_opt_alpha = cell(length(alpha), length(beta), length(gamma));
for i1 = 1:length(alpha)
    for i2 = 1:length(beta)
        for i3 = 1:length(gamma)
            input_dist_alpha{i1,i2,i3} = [(1-alpha(i1))*(1-beta(i2)) (1-alpha(i1))*beta(i2) alpha(i1)*gamma(i3) alpha(i1)*(1-gamma(i3))];
            rho = input_dist_alpha{i1,i2,i3}(1) + input_dist_alpha{i1,i2,i3}(4) - (input_dist_alpha{i1,i2,i3}(2) + input_dist_alpha{i1,i2,i3}(3));
            rho_mat = [1 rho; rho 1];
            
            sum_rate = @(w) 0;
            for i_user = 1:R
                W_mat = @(w) [w(1) w(2); w(3) w(4)];
                mean_vec = @(w) sqrt(P/trace(rho_mat*(W_mat(w)'*W_mat(w)))) * H(i_user,:) * W_mat(w) * s_all;
                sigma_vec = ones(1,2^R);
                alpha_vec = input_dist_alpha{i1,i2,i3};
                fun = entropy_function(mean_vec, sigma_vec, alpha_vec);
                entropy_Y = @(w) integral(@(y) fun(y,w), low_val, upper_val);

                idx_plus = (s_all(i_user,:) == 1);
                mean_vec_plus = @(w) sqrt(P/trace(rho_mat*(W_mat(w)'*W_mat(w)))) * H(i_user,:) * W_mat(w) * s_all(:,idx_plus);
                sigma_vec_plus = ones(1,sum(idx_plus));
                alpha_vec_plus = input_dist_alpha{i1,i2,i3}(idx_plus)/sum(input_dist_alpha{i1,i2,i3}(idx_plus));
                fun_plus = entropy_function(mean_vec_plus, sigma_vec_plus, alpha_vec_plus);
                idx_minus = (s_all(i_user,:) == -1);
                mean_vec_minus = @(w) sqrt(P/trace(rho_mat*(W_mat(w)'*W_mat(w)))) * H(i_user,:) * W_mat(w) * s_all(:,idx_minus);
                sigma_vec_minus = ones(1,sum(idx_minus));
                alpha_vec_minus = input_dist_alpha{i1,i2,i3}(idx_minus)/sum(input_dist_alpha{i1,i2,i3}(idx_minus));
                fun_minus = entropy_function(mean_vec_minus, sigma_vec_minus, alpha_vec_minus);
                entropy_YgivenX = @(w) sum(input_dist_alpha{i1,i2,i3}(idx_plus)) * integral(@(y) fun_plus(y,w), low_val, upper_val) + ...
                    sum(input_dist_alpha{i1,i2,i3}(idx_minus)) * integral(@(y) fun_minus(y,w), low_val, upper_val);

                sum_rate = @(w) sum_rate(w) + entropy_Y(w) - entropy_YgivenX(w);

                if i_user > 1
                    sum_rate = @(w) sum_rate(w) - h2(sum(input_dist_alpha{i1,i2,i3}(idx_plus)));
                end
            end
            
            idx_user1_plus = (s_all(1,:) == 1);
            idx_user1_minus = (s_all(1,:) == -1);
            input_dist_plus_givenX1 = input_dist_alpha{i1,i2,i3}(idx_user1_plus)/sum(input_dist_alpha{i1,i2,i3}(idx_user1_plus));
            input_dist_minus_givenX1 = input_dist_alpha{i1,i2,i3}(idx_user1_minus)/sum(input_dist_alpha{i1,i2,i3}(idx_user1_minus));
            sum_rate = @(w) sum_rate(w) + sum(input_dist_alpha{i1,i2,i3}(idx_user1_plus))*sum(-input_dist_plus_givenX1.*log2(input_dist_plus_givenX1)) + ...
                sum(input_dist_alpha{i1,i2,i3}(idx_user1_minus))*sum(-input_dist_minus_givenX1.*log2(input_dist_minus_givenX1));
            
            options = optimset('MaxFunEvals', 2000, 'MaxIter', 2000);
            w0 = W_initial_tr(:);
            [W_opt_alpha{i1,i2,i3}, opt_val] = fminsearch(@(w) -1*sum_rate(w), w0, options);
            rate_marton_alpha(i1,i2,i3) = -opt_val;
            
            if abs(alpha(i1)-0.5) < 1e-5 && abs(beta(i2)-0.5) < 1e-5 && abs(gamma(i3)-0.5) < 1e-5
                W_opt_symmetric = transpose(reshape(W_opt_alpha{i1,i2,i3}, R, []));
                sum_rate_symmetric = rate_marton_alpha(i1,i2,i3);
            end
        end
    end
end

[sum_rate, idx] = max(rate_marton_alpha, [], 'all', 'linear');
[i1_max, i2_max, i3_max] = ind2sub(size(rate_marton_alpha), idx);
input_dist = input_dist_alpha{i1_max, i2_max, i3_max};
W_opt = transpose(reshape(W_opt_alpha{i1_max, i2_max, i3_max}, R, []));
toc

end



function z = pdf_channel_output(y, w, mean_vec, sigma_vec, alpha_vec)
    z = 0 ;
    temp = mean_vec(w);
    for k = 1:length(sigma_vec)
        z = z + alpha_vec(k) .* 1/sqrt(2*pi*sigma_vec(k)^2) .* exp( -(y-temp(k)).^2/(2*sigma_vec(k)^2) );
    end
end

function z = entropy_function(mean_vec, sigma_vec, alpha_vec)
    z = @(y,w) max(0, -1 * pdf_channel_output(y, w, mean_vec, sigma_vec, alpha_vec) .* log2(pdf_channel_output(y, w, mean_vec, sigma_vec, alpha_vec)));
end
