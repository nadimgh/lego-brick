%testing script for the coding scheme of the uplink C-RAN problem

clearvars;
close all;

current_folder = pwd;
addpath(genpath([fileparts(current_folder), filesep, 'src' ]));

% If C<1, the input parameters should satisfy that (rate of compression + constraint_backoff) is smaller than C
% If C>=1, the input parameters should satisfy that rate of compression is smaller than C
% (if C>=1, we don't need to use constraint_backoff since there is no need for compression (recall the channel model includes a binary quantizer))

n = 1024;
sum_rate = 4/16;
rate_backoff = 4/32;        % rate backoff parameter for codes used in the channel coding part
constraint_backoff = 5/32;  % rate backoff parameter for codes used in the compression part
g = 0.9;
H = [1 g; g 1];
P_dB = -10:0.5:15;
C_vec = [1 1];
scl_flag = 0;
list_size = 1;
num_realizations = 5e4;
target_snr_flag = 1;

[rate_cell, dist_struct] = getMaximumAchievableUserRatesUplinkCRAN(H, P_dB, C_vec);
rate = getUserRatesUplinkCRAN(sum_rate, rate_backoff, constraint_backoff, rate_cell);

P = 10.^(P_dB/10);
[T, R] = size(H); % R-user, T-relay
[~, idx_target_snr] = min(abs(rate_cell{1}.max_sum_rate - sum_rate - R*rate_backoff));

checkInputParameters_UplinkCran(n, rate, dist_struct, H, C_vec+(C_vec(1)<1)*constraint_backoff);
k = n*rate;

matfile_capacity = sprintf('AchievableRates_UplinkCRAN_g%.2f_C%.2f.mat', g, C_vec(1));
if exist(matfile_capacity, 'file')
    res_capacity = load(matfile_capacity);
else
    error('Capacity matfile does not exist.');
end

nb_block_err = zeros(1, length(P));
nb_bit_err = zeros(1, length(P));
nb_block_err_user = zeros(R, length(P));
nb_bit_err_user = zeros(R, length(P));
nb_bit_err_relay = zeros(T, length(P));
avg_distortion_relay = zeros(T, length(P));
Pe_block = zeros(1, length(P));
Pe_bit = zeros(1, length(P));
Pe_block_user = zeros(R, length(P));
Pe_bit_user = zeros(R, length(P));
true_dist_input = zeros(length(P), 2^R);
true_dist_output_quantization = zeros(length(P), 2^(2*T));

tic
for i_p = 1:length(P)
    H_ch_effective = sqrt(P(i_p))*H;
    if target_snr_flag == 0
        asym_encoder1_dist = dist_struct.asym_encoder1(i_p,:);
        asym_encoder2_dist = dist_struct.asym_encoder2(i_p,:);
        asym_decoder1_dist = dist_struct.asym_decoder1(i_p,:);
        asym_decoder2_dist = dist_struct.asym_decoder2(i_p,:);
        lossy_encoder_dist = dist_struct.lossy_encoder(i_p,:);
        lossy_decoder_dist = dist_struct.lossy_decoder(i_p,:);
        wz_encoder_dist = dist_struct.wz_encoder(i_p,:);
        wz_decoder_dist = dist_struct.wz_decoder(i_p,:);
        partition = zeros(2,1); partition(1) = res_capacity.partition{1}(i_p); partition(2) = res_capacity.partition{2}(i_p);
    else
        asym_encoder1_dist = dist_struct.asym_encoder1(idx_target_snr,:);
        asym_encoder2_dist = dist_struct.asym_encoder2(idx_target_snr,:);
        asym_decoder1_dist = dist_struct.asym_decoder1(idx_target_snr,:);
        asym_decoder2_dist = dist_struct.asym_decoder2(idx_target_snr,:);
        lossy_encoder_dist = dist_struct.lossy_encoder(idx_target_snr,:);
        lossy_decoder_dist = dist_struct.lossy_decoder(idx_target_snr,:);
        wz_encoder_dist = dist_struct.wz_encoder(idx_target_snr,:);
        wz_decoder_dist = dist_struct.wz_decoder(idx_target_snr,:);
        partition = zeros(2,1); partition(1) = res_capacity.partition{1}(idx_target_snr); partition(2) = res_capacity.partition{2}(idx_target_snr);
    end
    
    bhatta_pcm = cell(R+T,2);
    bhatta_pcm{1,1} = computeBhattacharyyaPolar(n, asym_decoder1_dist);
    bhatta_pcm{1,2} = computeBhattacharyyaPolar(n, asym_encoder1_dist);
    bhatta_pcm{2,1} = computeBhattacharyyaPolar(n, asym_decoder2_dist);
    bhatta_pcm{2,2} = computeBhattacharyyaPolar(n, asym_encoder2_dist);
    bhatta_pcm{3,1} = computeBhattacharyyaPolar(n, lossy_encoder_dist);
    bhatta_pcm{3,2} = computeBhattacharyyaPolar(n, lossy_decoder_dist);
    bhatta_pcm{4,1} = computeBhattacharyyaPolar(n, wz_encoder_dist);
    bhatta_pcm{4,2} = computeBhattacharyyaPolar(n, wz_decoder_dist);
    
    H_struct = getParityCheckStruct(n, rate, bhatta_pcm);
    H_struct = StructToDouble(H_struct);
    
    for realization = 1:num_realizations
        if mod(realization, ceil(num_realizations/10)) == 1
            disp(['Sim iteration running = ', num2str(realization)]);
        end
        
        % Encoder at users
        u = cell(1,R);
        for i_user = 1:R
            u{i_user} = randi([0 1], k(i_user,1)-k(i_user,2), 1, 'int8');
        end
        [v_user, x_user] = UplinkCranEncoderUser(u, asym_encoder1_dist, asym_encoder2_dist, H_struct, 0, list_size);
        
        % Compute true input distribution
        x_mat = zeros(n, R);
        for i_user = 1:R
            x_mat(:,i_user) = x_user{i_user};
        end
        x_dec = bi2de(x_mat, 'left-msb')+1;
        [counts, indices] = groupcounts(x_dec);
        true_dist_input(i_p,indices) = true_dist_input(i_p,indices) + counts'/n;
        
        % Channel
        y_tilde = cell(1,T);
        for i_relay = 1:T
            y_tilde{i_relay} = zeros(n, 1);
            for i_user = 1:R
                y_tilde{i_relay} = y_tilde{i_relay} + H_ch_effective(i_relay, i_user) * (1-2*double(x_user{i_user}));
            end
            y_tilde{i_relay} = y_tilde{i_relay} + randn(n,1);
        end
        
        % Quantization
        y = cell(1,T);
        for i_relay = 1:T
            y{i_relay} = zeros(n, 1);
            y{i_relay}(y_tilde{i_relay} <= partition(i_relay)) = 1;
        end
        
        % Encoder at relays
        [v_relay, index_relay] = UplinkCranEncoderRelay(y, lossy_encoder_dist, wz_encoder_dist, H_struct, 0, list_size);
        
        % Decoder at CP
        [uhat, ~, yhat] = UplinkCranDecoderCP(index_relay, v_user, v_relay, lossy_decoder_dist, wz_decoder_dist, asym_decoder1_dist, asym_decoder2_dist, H_struct, scl_flag, list_size);
        
        % Compute true joint output-quantization distribution
        y_mat = zeros(n, 2*T);
        for i_relay = 1:T
            y_mat(:,i_relay) = y{i_relay};
        end
        for i_relay = T+1:2*T
            y_mat(:,i_relay) = yhat{i_relay-T};
        end
        y_dec = bi2de(y_mat, 'left-msb')+1;
        [counts, indices] = groupcounts(y_dec);
        true_dist_output_quantization(i_p,indices) = true_dist_output_quantization(i_p,indices) + counts'/n;
        
        % Error probability
        err = zeros(1, R);
        for i_user = 1:R
            err(i_user) = any(uhat{i_user} ~= u{i_user});
            if err(i_user)
                nb_block_err_user(i_user,i_p) = nb_block_err_user(i_user,i_p) + 1;
                nb_bit_err_user(i_user,i_p) = nb_bit_err_user(i_user,i_p) + sum(uhat{i_user} ~= u{i_user});
                nb_bit_err(i_p) = nb_bit_err(i_p) + sum(uhat{i_user} ~= u{i_user});
            end
        end
        
        for i_relay = 1:T
            nb_bit_err_relay(i_relay,i_p) = nb_bit_err_relay(i_relay,i_p) + sum(yhat{i_relay} ~= y{i_relay});
        end
        
        if any(err)
            nb_block_err(i_p) = nb_block_err(i_p) + 1;
        end
        
        if (nb_block_err(i_p) >= 100)
            break;
        end
    end
    
    Pe_block_user(:,i_p) = nb_block_err_user(:,i_p)/realization;
    Pe_block(i_p) = nb_block_err(i_p)/realization;
    avg_distortion_relay(:,i_p) = nb_bit_err_relay(:,i_p) / (realization * n);
    Pe_bit_user(:,i_p) = nb_bit_err_user(:,i_p) ./ (realization * (k(1:R,1) - k(1:R,2)));
    Pe_bit(i_p) = nb_bit_err(i_p) / (realization * sum(k(1:R,1) - k(1:R,2)));
    true_dist_input(i_p,:) = true_dist_input(i_p,:)/realization;
    true_dist_output_quantization(i_p,:) = true_dist_output_quantization(i_p,:)/realization;
    
    if (Pe_block(i_p) < 1e-4)
        break;
    end
end

toc

