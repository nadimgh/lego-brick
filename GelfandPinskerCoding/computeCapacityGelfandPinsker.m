function [beta_star, gamma_star, beta_star_bsc, capacity, capacity_bsc, capacity_symmetric] = computeCapacityGelfandPinsker(p_S, gauss_params)
%computes the capacity of a binary-input channel with state Y=X+S+Z, where
%S is binary and Z is Gaussian
%  p_S: 1x2 vector for the distribution of S (has the form [1-theta theta]) 
%  gauss_params: struct that includes (mean00, sigma00, alpha00, mean01, sigma01, alpha01, mean10, sigma10, alpha10, mean11, sigma11, alpha11) for the gaussian mixture parameters given s,x=0,1 (in this order)
%  beta_star: P(X=1|S=0) of the capacity-achieving conditional distribution
%  gamma_star: P(X=0|S=1) of the capacity-achieving conditional distribution
%  beta_star_bsc: crossover probability of the capacity-achieving BSC
%  capacity: capacity of the channel
%  capacity_bsc: capacity of the channel when p(x|s) is a BSC
%  capacity_symmetric: capacity of the channel when p(x) is Bern(0.5) and independent of S

mean00 = gauss_params.mean00; sigma00 = gauss_params.sigma00; alpha00 = gauss_params.alpha00;
mean01 = gauss_params.mean01; sigma01 = gauss_params.sigma01; alpha01 = gauss_params.alpha01;
mean10 = gauss_params.mean10; sigma10 = gauss_params.sigma10; alpha10 = gauss_params.alpha10;
mean11 = gauss_params.mean11; sigma11 = gauss_params.sigma11; alpha11 = gauss_params.alpha11;
mean_vec = [mean00 mean01 mean10 mean11];
sigma_vec = [sigma00 sigma01 sigma10 sigma11];
low_val = min(mean_vec)-5*max(sigma_vec);
upper_val = max(mean_vec)+5*max(sigma_vec);

theta = p_S(2);
beta = 0.01:0.01:0.5;
gamma = 0.01:0.01:0.5;

capacity_mat = zeros(length(beta), length(gamma));
capacity_symmetric_beta = zeros(1, length(beta));
for i = 1:length(beta)
    for j = 1:length(gamma)
        alpha_vec = [(1-theta)*(1-beta(i))*alpha00 (1-theta)*beta(i)*alpha10 theta*gamma(j)*alpha01 theta*(1-gamma(j))*alpha11];
        entropy_Y = integral(entropy_function(mean_vec, sigma_vec, alpha_vec), low_val, upper_val);
        
        % joint_dist_XS = [P(S=0,X=0) P(S=0,X=1) P(S=1,X=0) P(S=1,X=1)]
        joint_dist_XS = [(1-theta)*(1-beta(i)) (1-theta)*beta(i) theta*gamma(j) theta*(1-gamma(j))];
        p_X = transpose(sum(reshape(joint_dist_XS, 2, []), 2));
        p_SgivenX = joint_dist_XS./kron([1 1], p_X);
        entropy_YgivenX = p_X(1)*integral(entropy_function([mean00 mean10], [sigma00 sigma10], [p_SgivenX(1)*alpha00 p_SgivenX(3)*alpha10]), low_val, upper_val) + ...
            p_X(2)*integral(entropy_function([mean01 mean11], [sigma01 sigma11], [p_SgivenX(2)*alpha01 p_SgivenX(4)*alpha11]), low_val, upper_val);
        
        mutual_info_XS = computeMutualInformation(joint_dist_XS);
        capacity_mat(i,j) = entropy_Y - entropy_YgivenX - mutual_info_XS;
        
        if (beta(i) == gamma(j))
            capacity_symmetric_beta(i) = capacity_mat(i,j);
        end
    end
end
[capacity_beta, idx_beta] = max(capacity_mat);
[capacity, idx] = max(capacity_beta);
beta_star = beta(idx_beta(idx));
gamma_star = gamma(idx);

[capacity_bsc, idx] = max(capacity_symmetric_beta);
beta_star_bsc = beta(idx);

idx1 = find(beta == 0.5, 1);
idx2 = find(gamma == 0.5, 1);
capacity_symmetric = capacity_mat(idx1,idx2);

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
