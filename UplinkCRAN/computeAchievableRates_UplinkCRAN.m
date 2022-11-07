%compute the maximum achievable rates for the uplink C-RAN problem

clearvars;
close all;

current_folder = pwd;
addpath(genpath([fileparts(current_folder), filesep, 'src' ]));

g = 0.9;
H = [1 g; g 1];
P_dB = -10:0.5:15;
C_vec = [1 1];
[T, R] = size(H); % R-user, T-relay
P = 10.^(P_dB/10);

rng('shuffle');
tol = 1e-5;
s_all = 1-2*transpose(get_tuples(R));

% params = [P(X1=1) P(X2=1) P(Y1hat=1|Y1=0) P(Y1hat=1|Y1=1) P(Y2hat=1|Y2=0) P(Y2hat=1|Y2=1)]
input_dist = @(params) setInputDistFromParamsUplinkCRAN(params);
W_mat = @(power) [sqrt(power(1)) 0; 0 sqrt(power(2))];
x_all = @(power) W_mat(power) * s_all;
gauss_params = cell(1,T);
partition_fun = cell(1,T);
quantization_conditional_dist_user = cell(1,T); % i-th cell has [P(Yhat_i=1|X1=0,X2=0) P(Yhat_i=1|X1=0,X2=1) P(Yhat_i=1|X1=1,X2=0) P(Yhat_i=1|X1=1,X2=1)]
output_conditional_dist_user = cell(1,T); % i-th cell has [P(Y_i=1|X1=0,X2=0) P(Y_i=1|X1=0,X2=1) P(Y_i=1|X1=1,X2=0) P(Y_i=1|X1=1,X2=1)]
for i_relay = 1:T
    mean_vec = @(power) H(i_relay,:) * x_all(power);
    sigma_vec = ones(1,2^R);
    alpha_vec = input_dist;
    gauss_params{i_relay} = @(params,power) struct('alpha_vec', alpha_vec(params), 'mean_vec', mean_vec(power), 'sigma_vec', sigma_vec);
    partition_fun{i_relay} = @(params,power) getLloydMaxQuantizationParameters(gauss_params{i_relay}(params,power)); % Lloyd-Max algorithm to get optimal partition
    output_conditional_dist_user{i_relay} = @(params,power) 1 - Q((partition_fun{i_relay}(params,power)-mean_vec(power))./sigma_vec);
    quantization_conditional_dist_user{i_relay} = @(params,power) params(2*i_relay+1)*(1-output_conditional_dist_user{i_relay}(params,power)) ...
        + (1-params(2*i_relay+2))*output_conditional_dist_user{i_relay}(params,power);
end

output_dist_fun = @(params,power) setOutputDistFromConditionalDistUplinkCRAN(input_dist, output_conditional_dist_user, params, power); % returns [P(Y1=0,Y2=0) P(Y1=0,Y2=1) P(Y1=1,Y2=0) P(Y1=1,Y2=1)]
joint_output_quantization_dist_fun = @(params,power) setJointOutputQuantizationDistUplinkCRAN(output_dist_fun, params, power);
joint_input_quantization_dist_fun = @(params,power) setQuantizationDistFromConditionalDistUplinkCRAN(input_dist, quantization_conditional_dist_user, params, power);  % returns [P(Y1hat=0,Y2hat=0) P(Y1hat=0,Y2hat=1) P(Y1hat=1,Y2hat=0) P(Y1hat=1,Y2hat=1)]
quantization_dist_fun = @(params,power) transpose(sum(reshape(joint_input_quantization_dist_fun(params,power), 4, []), 2)); % returns [P(Y1hat=0,Y2hat=0) P(Y1hat=0,Y2hat=1) P(Y1hat=1,Y2hat=0) P(Y1hat=1,Y2hat=1)]

sum_rate_fun = @(params,power) computeEntropy(input_dist(params)) + computeEntropy(quantization_dist_fun(params,power)) - computeEntropy(joint_input_quantization_dist_fun(params,power));

tic
params = zeros(length(P), R+2*T);
input_dist = zeros(length(P), 2^R);
sum_rate = zeros(1, length(P));
joint_output_quantization_dist = zeros(length(P), 2^(2*T));
joint_input_quantization_dist = zeros(length(P), 2^(R+T));
partition = cell(1,T);
codebook = cell(1,T);
for i_relay = 1:T
    partition{i_relay} = zeros(length(P), 1);
    codebook{i_relay} = zeros(length(P), 2);
end

for i_p = 1:length(P)
    disp(['Sim iteration running = ', num2str(i_p)]);
    
    fun_obj = @(x) -1*sum_rate_fun(x(1:6), [P(i_p), P(i_p)]);
    [x_opt, fval] = patternsearch(fun_obj, [0.5 0.5 0.5*rand(1,4)], [], [], [], [], tol*ones(1,6), 0.5*ones(1,6), @(x)computeNonLinearConstraintFunction(x, H, P(i_p), C_vec));
    
    sum_rate(i_p) = -fval;
    params(i_p,:) = x_opt(1:6);
    input_dist(i_p,:) = setInputDistFromParamsUplinkCRAN(params(i_p,:));
    joint_output_quantization_dist(i_p,:) = joint_output_quantization_dist_fun(params(i_p,:), [P(i_p), P(i_p)]);
    joint_input_quantization_dist(i_p,:) = joint_input_quantization_dist_fun(params(i_p,:), [P(i_p), P(i_p)]);
    for i_relay = 1:T
        [partition{i_relay}(i_p), codebook{i_relay}(i_p,:), ~] = getLloydMaxQuantizationParameters(gauss_params{i_relay}(params(i_p,:), [P(i_p), P(i_p)]));
    end
end

matfile = sprintf('Results_Capacity/AchievableRates_UplinkCRAN_g%.2f_C%.2f.mat', g, C_vec(1));
if ~exist(matfile, 'file')
    save(matfile, 'params', 'input_dist', 'sum_rate', 'P_dB', 'H', 'C_vec', 'joint_output_quantization_dist', 'joint_input_quantization_dist', 'partition', 'codebook');
end
toc
