function [input_dist, sum_rate, sum_rate_symmetric] = computeCapacityMartonWithPrecoding(H, P, W_precoding)
%compute the achievable sum rate of marton coding when followed by
%precoding using the matrix W_precoding

[R, ~] = size(H);
alpha = 0.01:0.01:0.5;
beta = 0.01:0.01:0.5;
gamma = 0.01:0.01:0.5;

s_all = 1-2*transpose(get_tuples(R));
x_all = (sqrt(P)/norm(W_precoding(:))) * W_precoding * s_all; % E[xx^T] = I in this case, so average transmitted power is norm(W_precoding(:))^2

low_val = min(H * x_all, [], 'all')-5;
upper_val = max(H * x_all, [], 'all')+5;

input_dist_alpha = cell(length(alpha), length(beta), length(gamma));
rate_marton_alpha = zeros(length(alpha), length(beta), length(gamma));
for i1 = 1:length(alpha)
    for i2 = 1:length(beta)
        for i3 = 1:length(gamma)
            input_dist_alpha{i1,i2,i3} = [(1-alpha(i1))*(1-beta(i2)) (1-alpha(i1))*beta(i2) alpha(i1)*gamma(i3) alpha(i1)*(1-gamma(i3))];
            for i_user = 1:R
                mean_vec = H(i_user,:)*x_all;
                sigma_vec = ones(1,2^R);
                alpha_vec = input_dist_alpha{i1,i2,i3};
                entropy_Y = integral(entropy_function(mean_vec, sigma_vec, alpha_vec), low_val, upper_val);

                idx_plus = (s_all(i_user,:) == 1);
                mean_vec_plus = H(i_user,:)*x_all(:,idx_plus);
                sigma_vec_plus = ones(1,length(idx_plus));
                alpha_vec_plus = input_dist_alpha{i1,i2,i3}(idx_plus)/sum(input_dist_alpha{i1,i2,i3}(idx_plus));
                idx_minus = (s_all(i_user,:) == -1);
                mean_vec_minus = H(i_user,:)*x_all(:,idx_minus);
                sigma_vec_minus = ones(1,length(idx_minus));
                alpha_vec_minus = input_dist_alpha{i1,i2,i3}(idx_minus)/sum(input_dist_alpha{i1,i2,i3}(idx_minus));
                entropy_YgivenX = sum(input_dist_alpha{i1,i2,i3}(idx_plus)) * integral(entropy_function(mean_vec_plus, sigma_vec_plus, alpha_vec_plus), low_val, upper_val) + ...
                    sum(input_dist_alpha{i1,i2,i3}(idx_minus)) * integral(entropy_function(mean_vec_minus, sigma_vec_minus, alpha_vec_minus), low_val, upper_val);

                rate_marton_alpha(i1,i2,i3) = rate_marton_alpha(i1,i2,i3) + entropy_Y - entropy_YgivenX;

                if i_user > 1
                    rate_marton_alpha(i1,i2,i3) = rate_marton_alpha(i1,i2,i3) - h2(sum(input_dist_alpha{i1,i2,i3}(idx_plus)));
                end
            end
            
            idx_user1_plus = (s_all(1,:) == 1);
            idx_user1_minus = (s_all(1,:) == -1);
            input_dist_plus_givenX1 = input_dist_alpha{i1,i2,i3}(idx_user1_plus)/sum(input_dist_alpha{i1,i2,i3}(idx_user1_plus));
            input_dist_minus_givenX1 = input_dist_alpha{i1,i2,i3}(idx_user1_minus)/sum(input_dist_alpha{i1,i2,i3}(idx_user1_minus));
            rate_marton_alpha(i1,i2,i3) = rate_marton_alpha(i1,i2,i3) + sum(input_dist_alpha{i1,i2,i3}(idx_user1_plus))*sum(-input_dist_plus_givenX1.*log2(input_dist_plus_givenX1)) + ...
                sum(input_dist_alpha{i1,i2,i3}(idx_user1_minus))*sum(-input_dist_minus_givenX1.*log2(input_dist_minus_givenX1));
            
            if abs(alpha(i1)-0.5) < 1e-5 && abs(beta(i2)-0.5) < 1e-5 && abs(gamma(i3)-0.5) < 1e-5
                sum_rate_symmetric = rate_marton_alpha(i1,i2,i3);
            end
        end
    end
end

[sum_rate, idx] = max(rate_marton_alpha, [], 'all', 'linear');
[i1_max, i2_max, i3_max] = ind2sub(size(rate_marton_alpha), idx);
input_dist = input_dist_alpha{i1_max, i2_max, i3_max};

end


function z = pdf_channel_output(y, mean, sigma, alpha)
    z = 0 ;
    for k = 1:length(mean)
        z = z + alpha(k) .* 1/sqrt(2*pi*sigma(k)^2) .* exp( -(y-mean(k)).^2/(2*sigma(k)^2) );
    end
end

function z = entropy_function(mean, sigma, alpha)
    z = @(y) -1 * pdf_channel_output(y, mean, sigma, alpha) .* log2(pdf_channel_output(y, mean, sigma, alpha));
end

