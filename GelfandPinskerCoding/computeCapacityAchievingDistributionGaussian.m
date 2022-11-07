function [input_dist, input_dist_bsc, capacity, capacity_bsc, capacity_symmetric] = computeCapacityAchievingDistributionGaussian(p_s, g, P_dB)
%compute the capacity-achieving distributions for the Gaussian
%channel Y=X+gS+Z, such that (X,S) are each subject to a power constraint
%P, where
%  p_s: distribution of S (has the form [1-theta theta])
%  g: channel parameter
%  P_dB: power constraint in dB
%  save_flag: save output to a mat file
%  input_dist: capacity-achieving input distribution p(s,x) (in this order)
%  input_dist_bsc: capacity-achieving input distribution when p(x|s) is a bsc
%  capacity: channel capacity
%  capacity_bsc: channel capacity when p(x|s) is a bsc
%  capacity_symmetric: channel capacity when p(x) is Bern(0.5) and independent of S

P = 10.^(P_dB/10);
theta = p_s(2);
alpha00 = 1; alpha01 = 1; alpha10 = 1; alpha11 = 1;
sigma00 = 1; sigma01 = 1; sigma10 = 1; sigma11 = 1;

capacity = zeros(1, length(P));
capacity_bsc = zeros(1, length(P));
capacity_symmetric = zeros(1, length(P));
beta_star = zeros(length(P), 1);
gamma_star = zeros(length(P), 1);
beta_star_bsc = zeros(length(P), 1);
for i = 1:length(P)
    mean00 = sqrt(P(i))*(1+g); mean01 = sqrt(P(i))*(-1+g); mean10 = sqrt(P(i))*(1-g); mean11 = sqrt(P(i))*(-1-g);
    gauss_params = struct('mean00', mean00, 'mean01', mean01, 'mean10', mean10, 'mean11', mean11, 'sigma00', sigma00, 'sigma01', sigma01, 'sigma10', sigma10, 'sigma11', sigma11, ...
        'alpha00', alpha00, 'alpha01', alpha01, 'alpha10', alpha10, 'alpha11', alpha11);
    [beta_star(i), gamma_star(i), beta_star_bsc(i), capacity(i), capacity_bsc(i), capacity_symmetric(i)] = computeCapacityGelfandPinsker(p_s, gauss_params);
end

input_dist = [(1-theta)*(1-beta_star) (1-theta)*beta_star theta*gamma_star theta*(1-gamma_star)];
input_dist_bsc = [(1-theta)*(1-beta_star_bsc) (1-theta)*beta_star_bsc theta*beta_star_bsc theta*(1-beta_star_bsc)];

end

