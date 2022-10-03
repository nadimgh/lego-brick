function [sum_rate, rate_user] = computeCapacityPrecoding(H, P, W_precoding)
%compute the achievable sum rate of applying precoding matrix W_precoding

[R, ~] = size(H);

s_all = 1-2*transpose(get_tuples(R));
x_all = (sqrt(P)/norm(W_precoding(:))) * W_precoding * s_all;

low_val = min(H * x_all, [], 'all')-5;
upper_val = max(H * x_all, [], 'all')+5;
rate_user = zeros(1, R);

for i_user = 1:R
    mean_vec = H(i_user,:)*x_all;
    sigma_vec = ones(1,2^R);
    alpha_vec = (1/(2^R))*ones(1,2^R);
    entropy_Y = integral(entropy_function(mean_vec, sigma_vec, alpha_vec), low_val, upper_val);

    idx_plus = (s_all(i_user,:) == 1);
    mean_vec_plus = H(i_user,:)*x_all(:,idx_plus);
    sigma_vec_plus = ones(1,length(idx_plus));
    alpha_vec_plus = (1/(2^(R-1)))*ones(1,length(idx_plus));
    idx_minus = (s_all(i_user,:) == -1);
    mean_vec_minus = H(i_user,:)*x_all(:,idx_minus);
    sigma_vec_minus = ones(1,length(idx_minus));
    alpha_vec_minus = (1/(2^(R-1)))*ones(1,length(idx_minus));
    entropy_YgivenX = 0.5 * integral(entropy_function(mean_vec_plus, sigma_vec_plus, alpha_vec_plus), low_val, upper_val) + ...
        0.5 * integral(entropy_function(mean_vec_minus, sigma_vec_minus, alpha_vec_minus), low_val, upper_val);

    rate_user(i_user) = entropy_Y - entropy_YgivenX;
end

sum_rate = sum(rate_user);

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

