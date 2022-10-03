%testing script for lossy source coding scheme for an asymmetric source

clearvars;
close all;

current_folder = pwd;
addpath(genpath([fileparts(current_folder), filesep, 'src' ]));

n = 1024;
rate_backoff = 1/8;
theta = 0.3;
distortion_level = 0.02;%:0.02:0.28;
scl_flag = 0;
list_size = 1;
num_realizations = 1e3;

rate = zeros(length(distortion_level), 2);
nb_bit_err = zeros(1, length(distortion_level));
nb_block_err_lossless = zeros(1, length(distortion_level));
nb_bit_err_lossless = zeros(1, length(distortion_level));
alpha_encoder = zeros(1, length(distortion_level));
nb_bit_err_encoder = zeros(1, length(distortion_level));
alpha_achieved = zeros(1, length(distortion_level));

tic
for i_d = 1:length(distortion_level)
    alpha = (theta-distortion_level(i_d))/(1-2*distortion_level(i_d));
    target_dist = [(1-alpha)*(1-distortion_level(i_d)), alpha*distortion_level(i_d), (1-alpha)*distortion_level(i_d), alpha*(1-distortion_level(i_d))];
    target_dist = target_dist(:);
    p_xhat_target = [1-alpha alpha];
    
    bhatta_pcm = cell(1,2);
    bhatta_pcm{1} = computeBhattacharyyaPolar(n, target_dist);
    bhatta_pcm{2} = computeBhattacharyyaPolar(n, p_xhat_target);
    
    rate(i_d,:) = getLossySourceCodeRates(target_dist, rate_backoff);
    if any(abs(n*rate(i_d,:) - floor(n*rate(i_d,:))) > 1e-5, 'all')
        error('n*rate should be an integer.');
    end
    
    H_struct = getParityCheckStruct(n, rate(i_d,:), bhatta_pcm);
    H_struct = StructToDouble(H_struct);
    H_struct = H_struct{1};
    
    for realization = 1:num_realizations
        if mod(realization, ceil(num_realizations/10)) == 1
            disp(['Sim iteration running = ', num2str(realization)]);
        end

        x = randsmpl([1-theta theta], n, 1, 'int8')-1;    % independent sampling from p_x 

        [index, v, x_encoder] = LossySourceEncoder(x, target_dist, H_struct, 0, list_size); % no need for SCL decoding for shaping
        xhat = LossySourceDecoder(index, v, p_xhat_target, H_struct, scl_flag, list_size);

        nb_bit_err(i_d) = nb_bit_err(i_d) + sum(xhat ~= x);
        alpha_encoder(i_d) = alpha_encoder(i_d) + sum(x_encoder)/n;
        nb_bit_err_encoder(i_d) = nb_bit_err_encoder(i_d) + sum(x ~= x_encoder);
        alpha_achieved(i_d) = alpha_achieved(i_d) + sum(xhat)/n;
        nb_block_err_lossless(i_d) = nb_block_err_lossless(i_d) + any(xhat ~= x_encoder);
        nb_bit_err_lossless(i_d) = nb_bit_err_lossless(i_d) + sum(xhat ~= x_encoder);
    end
end

avg_distortion = nb_bit_err/(realization*n);
alpha_encoder = alpha_encoder/realization;
avg_distortion_encoder = nb_bit_err_encoder/(realization*n);
alpha_achieved = alpha_achieved/realization;
Pe_block_lossless = nb_block_err_lossless/realization;
Pe_bit_lossless = nb_bit_err_lossless/(realization*n);
true_dist_encoder = [transpose((1-alpha_encoder).*(1-avg_distortion_encoder)), transpose(alpha_encoder.*avg_distortion_encoder), transpose((1-alpha_encoder).*avg_distortion_encoder), transpose(alpha_encoder.*(1-avg_distortion_encoder))];
true_dist = [transpose((1-alpha_achieved).*(1-avg_distortion)), transpose(alpha_achieved.*avg_distortion), transpose((1-alpha_achieved).*avg_distortion), transpose(alpha_achieved.*(1-avg_distortion))];

toc

