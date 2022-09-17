function [sum_rate, user_rate] = computeCapacityMartonSpecificInputDist(H, P, W_precoding, input_dist)
%compute the maximum achievable sum-rate (overall and per user) for a specific input distribution and
%precoding matrix

[R, ~] = size(H);

s_all = 1-2*transpose(get_tuples(R));
rho = input_dist(1) + input_dist(4) - (input_dist(2) + input_dist(3));
rho_mat = [1 rho; rho 1];
x_all = sqrt(P/trace(rho_mat*(W_precoding'*W_precoding))) * W_precoding * s_all;

low_val = min(H * x_all, [], 'all')-5;
upper_val = max(H * x_all, [], 'all')+5;

user_rate = zeros(1, R);
for i_user = 1:R
    mean_vec = H(i_user,:)*x_all;
    sigma_vec = ones(1,2^R);
    alpha_vec = input_dist;
    entropy_Y = integral(entropy_function(mean_vec, sigma_vec, alpha_vec), low_val, upper_val);

    idx_plus = (s_all(i_user,:) == 1);
    mean_vec_plus = H(i_user,:)*x_all(:,idx_plus);
    sigma_vec_plus = ones(1,length(idx_plus));
    alpha_vec_plus = input_dist(idx_plus)/sum(input_dist(idx_plus));
    idx_minus = (s_all(i_user,:) == -1);
    mean_vec_minus = H(i_user,:)*x_all(:,idx_minus);
    sigma_vec_minus = ones(1,length(idx_minus));
    alpha_vec_minus = input_dist(idx_minus)/sum(input_dist(idx_minus));
    entropy_YgivenX = sum(input_dist(idx_plus)) * integral(entropy_function(mean_vec_plus, sigma_vec_plus, alpha_vec_plus), low_val, upper_val) + ...
        sum(input_dist(idx_minus)) * integral(entropy_function(mean_vec_minus, sigma_vec_minus, alpha_vec_minus), low_val, upper_val);

    if i_user == 1
        user_rate(i_user) = entropy_Y - entropy_YgivenX;
    elseif i_user < R
        input_dist_user = sum(reshape(input_dist, 2^(R-i_user), []));
        user_rate(i_user) = entropy_Y - entropy_YgivenX - computeMutualInformation(input_dist_user);
    else
        user_rate(i_user) = entropy_Y - entropy_YgivenX - computeMutualInformation(input_dist);
    end
end
            
sum_rate = sum(user_rate);

end


function z = pdf_channel_output(y, mean, sigma, alpha)
    z = 0 ;
    for k = 1:length(mean)
        z = z + alpha(k) .* 1/sqrt(2*pi*sigma(k)^2) .* exp( -(y-mean(k)).^2/(2*sigma(k)^2) );
    end
end

function z = entropy_function(mean, sigma, alpha)
    z = @(y) max(0, -1 * pdf_channel_output(y, mean, sigma, alpha) .* log2(pdf_channel_output(y, mean, sigma, alpha)));
end

