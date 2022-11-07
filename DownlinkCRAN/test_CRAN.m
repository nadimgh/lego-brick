%testing script for the coding scheme of the downlink C-RAN problem

clearvars;
close all;

current_folder = pwd;
addpath(genpath([fileparts(current_folder), filesep, 'src' ]));

% If C<1, the input parameters should satisfy that (rate of compression + constraint_backoff) is smaller than C
% If C>=1, the input parameters should satisfy that rate of compression is smaller than C 
% (if C>=1, we don't need to use constraint_backoff since there is no need for compression)

n = 1024;
sum_rate = 12/16;
rate_backoff = 4/32;        % rate backoff parameter for codes used in the channel coding part
constraint_backoff = 5/32;  % rate backoff parameter for codes used in the compression part
g = 0.9;
H = [1 g; g 1];
P_dB = -10:0.5:15;
C_vec = [1 1];
scl_flag = 0;
list_size = 1;
num_realizations = 5e4;
P = 10.^(P_dB/10);
[R, T] = size(H);
sum_power_flag = 0;
target_snr_flag = 1;

rate_cell = getMaximumAchievableUserRatesDownlinkCRAN(H, P_dB, 'cran', C_vec, sum_power_flag);
rate = getUserRatesDownlinkCRAN(sum_rate, rate_backoff, constraint_backoff, rate_cell, 'cran');
[~, idx_target_snr] = min(abs(rate_cell{1}.rate_cran - sum_rate - R*rate_backoff));

target_dist = rate_cell{1}.input_dist;
W_precoding = rate_cell{3};
k = n*rate;

checkInputParameters_DownlinkCRAN(n, rate, target_dist, H, C_vec+(C_vec<1)*constraint_backoff); 

x_all = 1-2*transpose(get_tuples(R+T));

nb_block_err = zeros(1, length(P));
nb_bit_err = zeros(1, length(P));
nb_block_err_user = zeros(R, length(P));
nb_bit_err_user = zeros(R, length(P));
nb_block_err_relay = zeros(T, length(P));
nb_bit_err_relay = zeros(T, length(P));
Pe_block = zeros(1, length(P));
Pe_bit = zeros(1, length(P));
Pe_block_user = zeros(R, length(P));
Pe_bit_user = zeros(R, length(P));
Pe_block_relay = zeros(T, length(P));
Pe_bit_relay = zeros(T, length(P));
true_dist_cp = zeros(length(P), 2^(R+T));

tic
for i_p = 1:length(P)
    H_ch_effective = H*W_precoding{i_p};
    if target_snr_flag == 0
        target_dist_encoder = target_dist(i_p,:);
    else
        target_dist_encoder = target_dist(idx_target_snr,:);
    end
    
    gauss_params = cell(1,R);
    bhatta_pcm = cell(R+T,2);
    p_x_target = cell(1,R+T); % i-th element is p(u_1^i) or p(u_1^R,x_1^(i-R))
    p_x_user_cp = cell(1,R); % i-th element is p(u_i)
    p_x_user_relay = cell(1,T); % i-th element is p(x_i)
    for i_user = 1:R
        idx0 = (x_all(i_user,:) == 1);
        mean0 = H_ch_effective(i_user,:) * x_all(R+1:R+T,idx0);
        alpha0 = target_dist_encoder(idx0)/sum(target_dist_encoder(idx0));
        sigma0 = ones(1, 2^(R+T-1));
        
        idx1 = (x_all(i_user,:) == -1);
        mean1 = H_ch_effective(i_user,:) * x_all(R+1:R+T,idx1);
        alpha1 = target_dist_encoder(idx1)/sum(target_dist_encoder(idx1));
        sigma1 = ones(1, 2^(R+T-1));
        
        p_x_user_cp{i_user} = [sum(target_dist_encoder(idx0)) sum(target_dist_encoder(idx1))];
        
        gauss_params{i_user} = struct('p_x', p_x_user_cp{i_user}, 'mean0', mean0, 'alpha0', alpha0, 'sigma0', sigma0, 'mean1', mean1, 'alpha1', alpha1, 'sigma1', sigma1);
        
        p_x_target{i_user} = sum(reshape(target_dist_encoder, 2^(R+T-i_user), []));
        
        bhatta_pcm{i_user,1} = computeBhattacharyyaPolarGaussian(n, p_x_user_cp{i_user}, gauss_params{i_user}, 1);
        bhatta_pcm{i_user,2} = computeBhattacharyyaPolar(n, p_x_target{i_user});
    end
    
    for i_user = R+1:R+T
        idx0 = (x_all(i_user,:) == 1);
        idx1 = (x_all(i_user,:) == -1);
        p_x_user_relay{i_user-R} = [sum(target_dist_encoder(idx0)) sum(target_dist_encoder(idx1))];
        if i_user < R+T
            p_x_target{i_user} = sum(reshape(target_dist_encoder, 2^(R+T-i_user), []));
        else
            p_x_target{i_user} = target_dist_encoder;
        end
        bhatta_pcm{i_user,1} = computeBhattacharyyaPolar(n, p_x_target{i_user});
        bhatta_pcm{i_user,2} = computeBhattacharyyaPolar(n, p_x_user_relay{i_user-R});
    end
    
    H_struct = getParityCheckStruct(n, rate, bhatta_pcm);
    H_struct = StructToDouble(H_struct);
    
    for realization = 1:num_realizations
        if mod(realization, ceil(num_realizations/10)) == 1
            disp(['Sim iteration running = ', num2str(realization)]);
        end
        
        % Encoder at CP
        u = cell(1,R);
        for i_user = 1:R
            u{i_user} = randi([0 1], k(i_user,1)-k(i_user,2), 1, 'int8');
        end
        [index, v, x_cp] = CranEncoderCP(u, target_dist_encoder, H_struct, 0, list_size);
        
        % Compute true input distribution at central processor
        x_mat = zeros(n, R+T);
        for i_user = 1:R+T
            x_mat(:,i_user) = x_cp{i_user};
        end
        x_dec = bi2de(x_mat, 'left-msb')+1;
        [counts, indices] = groupcounts(x_dec);
        true_dist_cp(i_p,indices) = true_dist_cp(i_p,indices) + counts'/n;
        
        % Encoder at relays
        err_relay = zeros(1,T);
        x_relay = CranEncoderRelay(index, v, p_x_user_relay, H_struct, scl_flag, list_size);
        for i_user = R+1:R+T
            err_relay(i_user-R) = any(x_relay{i_user-R} ~= x_cp{i_user});
            if err_relay(i_user-R)
                nb_block_err_relay(i_user-R,i_p) = nb_block_err_relay(i_user-R,i_p) + 1;
                nb_bit_err_relay(i_user-R,i_p) = nb_bit_err_relay(i_user-R,i_p) + sum(x_relay{i_user-R} ~= x_cp{i_user});
            end
        end
        
        % Channel
        y = cell(1,R);
        for i_user = 1:R
            y{i_user} = zeros(n, 1);
            for i_relay = 1:T
                y{i_user} = y{i_user} + H_ch_effective(i_user, i_relay) * (1-2*double(x_relay{i_relay}));
            end
            y{i_user} = y{i_user} + randn(n,1);
        end
        
        % Decoder
        uhat = CranDecoder(y, v, H_struct, gauss_params, scl_flag, list_size);
        
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
        
        if any(err)
            nb_block_err(i_p) = nb_block_err(i_p) + 1;
        end
        
        if (nb_block_err(i_p) >= 100)
            break;
        end
    end
    
    Pe_block_relay(:,i_p) = nb_block_err_relay(:,i_p)/realization;
    Pe_block_user(:,i_p) = nb_block_err_user(:,i_p)/realization;
    Pe_block(i_p) = nb_block_err(i_p)/realization;
    Pe_bit_relay(:,i_p) = nb_bit_err_relay(:,i_p) / (realization * n);
    Pe_bit_user(:,i_p) = nb_bit_err_user(:,i_p) ./ (realization * (k(1:R,1) - k(1:R,2)));
    Pe_bit(i_p) = nb_bit_err(i_p) / (realization * sum(k(1:R,1) - k(1:R,2)));
    true_dist_cp(i_p,:) = true_dist_cp(i_p,:)/realization;
    
    if (Pe_block(i_p) < 1e-4)
        break;
    end
end

toc

