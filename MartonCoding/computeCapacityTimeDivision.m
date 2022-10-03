function sum_rate = computeCapacityTimeDivision(H, P)
%compute the achievable sum rate of the time-division coding scheme, where
%only one user accesses the channel

[R, ~] = size(H);

P_opt_total = zeros(1, R); % total optimal power when data is sent to a specific user
P_opt_user = zeros(R); % optimal power at each antenna when data is sent to a specific user
for i_user = 1:R
    den = sum(H(i_user,:).^2);
    P_opt_user(i_user,:) = H(i_user,:).^2*P/den;
end
for i_user = 1:R
    P_opt_total(i_user) = (H(i_user,:)*transpose(sqrt(P_opt_user(i_user,:))))^2;
end
P_opt = max(P_opt_total);

low_val = min(sqrt(P_opt) * H * [1 -1; 1 -1], [], 'all')-5;
upper_val = max(sqrt(P_opt) * H * [1 -1; 1 -1], [], 'all')+5;

mean_vec = [sqrt(P_opt) -sqrt(P_opt)];
sigma_vec = ones(1,2);
alpha_vec = 0.5*ones(1,2);
entropy_Y = integral(entropy_function(mean_vec, sigma_vec, alpha_vec), low_val, upper_val);

mean_vec_plus = sqrt(P_opt);
mean_vec_minus = -sqrt(P_opt);
entropy_YgivenX = 0.5 * integral(entropy_function(mean_vec_plus, 1, 1), low_val, upper_val) + ...
    0.5 * integral(entropy_function(mean_vec_minus, 1, 1), low_val, upper_val);

sum_rate = entropy_Y - entropy_YgivenX;

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
