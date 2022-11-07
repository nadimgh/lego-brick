%testing script for a Gelfand-Pinsker coding scheme over the Gaussian channel Y=X+gS+Z, where X and S are binary and belong to {+/-sqrt(P)}

clearvars;
close all;

current_folder = pwd;
addpath(genpath([fileparts(current_folder), filesep, 'src' ]));

n = 1024;
overall_rate = 1/2;
rate_backoff = 3/16;
P_dB = -10:0.5:15;
g = 0.9;
theta = 0.1;
scl_flag = 0;
list_size = 1;
target_snr_flag = 1;
num_realizations = 5e4;

p_s = [1-theta theta];
[rate, target_dist, idx_target_snr] = getGelfandPinskerCodeRates(overall_rate, rate_backoff, p_s, g, P_dB);

k = round(n*rate); k1 = k(1); k2 = k(2);
P = 10.^(P_dB/10);

sigma00 = 1; sigma01 = 1; sigma10 = 1; sigma11 = 1;
alpha00 = 1; alpha01 = 1; alpha10 = 1; alpha11 = 1;

nb_block_err = zeros(1, length(P));
nb_bit_err = zeros(1, length(P));
Pe_block = zeros(1, length(P));
Pe_bit = zeros(1, length(P));
true_dist = zeros(length(P), 4);

tic
for i_p = 1:length(P)
    if target_snr_flag == 0
        target_dist_encoder = target_dist(i_p,:);
    else
        target_dist_encoder = target_dist(idx_target_snr,:);
    end
    p_x = transpose(sum(reshape(target_dist_encoder, 2, []), 2));
    beta_tilde = target_dist_encoder(2)/p_x(2);  % P(S=0 | X=1)
    gamma_tilde = target_dist_encoder(3)/p_x(1); % P(S=1 | X=0)
    
    mean00 = sqrt(P(i_p))*(1+g); mean01 = sqrt(P(i_p))*(-1+g); mean10 = sqrt(P(i_p))*(1-g); mean11 = sqrt(P(i_p))*(-1-g);
    gauss_params = struct('p_x', target_dist_encoder, 'mean00', mean00, 'mean01', mean01, 'mean10', mean10, 'mean11', mean11, 'sigma00', sigma00, 'sigma01', sigma01, 'sigma10', sigma10, 'sigma11', sigma11, ...
        'alpha00', alpha00, 'alpha01', alpha01, 'alpha10', alpha10, 'alpha11', alpha11);
    gauss_params_XY = struct('p_x', p_x, 'mean0', [mean00 mean10], 'sigma0', [sigma00 sigma10], 'alpha0', [(1-gamma_tilde)*alpha00 gamma_tilde*alpha10], ...
        'mean1', [mean01 mean11], 'sigma1', [sigma01 sigma11], 'alpha1', [beta_tilde*alpha01 (1-beta_tilde)*alpha11]);
    
    bhatta_pcm = cell(1,2);
    bhatta_pcm{1} = computeBhattacharyyaPolarGaussian(n, p_x, gauss_params_XY, 1);
    bhatta_pcm{2} = computeBhattacharyyaPolar(n, p_x);
    H_struct = getParityCheckStruct(n, rate, bhatta_pcm);
    H_struct = StructToDouble(H_struct);
    H_struct = H_struct{1};
    
    for realization = 1:num_realizations
        if mod(realization, ceil(num_realizations/10)) == 1
            disp(['Sim iteration running = ', num2str(realization)]);
        end

        u = randi([0 1], k1-k2, 1, 'int8');
        s = randsmpl(p_s, n, 1, 'int8')-1;    % independent sampling from p_s

        [x, v] = GelfandPinskerEncoder(u, s, target_dist_encoder, H_struct, 0, list_size);
        true_dist(i_p,1) = true_dist(i_p,1) + sum(s==0 & x==0)/n;
        true_dist(i_p,2) = true_dist(i_p,2) + sum(s==0 & x==1)/n;
        true_dist(i_p,3) = true_dist(i_p,3) + sum(s==1 & x==0)/n;
        true_dist(i_p,4) = true_dist(i_p,4) + sum(s==1 & x==1)/n;

        y = generateGaussianMixtureChannelOutput(x, s, gauss_params);

        uhat = GelfandPinskerDecoder(y, v, target_dist_encoder, H_struct, 1, gauss_params_XY, scl_flag, list_size);

        err = any(uhat ~= u);
        if err
            nb_block_err(i_p) = nb_block_err(i_p) + 1;
            nb_bit_err(i_p) = nb_bit_err(i_p) + sum(uhat ~= u);
        end

        if (nb_block_err(i_p) >= 100)
            break;
        end
    end
    
    Pe_block(i_p) = nb_block_err(i_p)/realization;
    Pe_bit(i_p) = nb_bit_err(i_p)/(realization*(k1-k2));
    true_dist(i_p,:) = true_dist(i_p,:)/realization;
    
    if (Pe_block(i_p) < 1e-4)
        break;
    end
end

toc

