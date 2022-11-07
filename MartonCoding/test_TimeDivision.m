%testing script for time-division over a two-user broadcast channel 

clearvars;
close all;

current_folder = pwd;
addpath(genpath([fileparts(current_folder), filesep, 'src' ]));

n = 1024;
sum_rate = 1;
g = 0.9;
H = [1 g; g 1];
P_dB = -10:0.5:15;
scl_flag = 0;
list_size = 1;
num_realizations = 5e4;

[R, ~] = size(H);
P = 10.^(P_dB/10);

rate_struct = getMaximumAchievableUserRates(H, P_dB, 'time_division');
rate = getUserRates(sum_rate, rate_struct, 'time_division');
k = n*rate;

P_opt_total = zeros(R, length(P)); % total optimal power when data is sent to a specific user
P_opt_user = cell(1, length(P)); % optimal power at each antenna when data is sent to a specific user
for i_p = 1:length(P)
    P_opt_user{i_p} = zeros(R);
    for i_user = 1:R
        den = sum(H(i_user,:).^2);
        P_opt_user{i_p}(i_user,:) = H(i_user,:).^2*P(i_p)/den;
    end
    for i_user = 1:R
        P_opt_total(i_user,i_p) = (H(i_user,:)*transpose(sqrt(P_opt_user{i_p}(i_user,:))))^2;
    end
end

nb_block_err = zeros(1, length(P));
nb_bit_err = zeros(1, length(P));
Pe_block = zeros(1, length(P));
Pe_bit = zeros(1, length(P));

tic
for i_p = 1:length(P)
    gauss_params = cell(1,R);
    info_set = cell(1,R);
    for i_user = 1:R
        mean0 = sqrt(P_opt_total(i_user,i_p));
        alpha0 = 1;
        sigma0 = 1;
        
        mean1 = -sqrt(P_opt_total(i_user,i_p));
        alpha1 = 1;
        sigma1 = 1;
        
        gauss_params{i_user} = struct('p_x', [0.5 0.5], 'mean0', mean0, 'alpha0', alpha0, 'sigma0', sigma0, 'mean1', mean1, 'alpha1', alpha1, 'sigma1', sigma1);
        
        bhatta = computeBhattacharyyaPolarGaussian(n, [0.5 0.5], gauss_params{i_user}, 0);
        [~, idxBhatta] = sort(bhatta);
        info_set{i_user} = zeros(1,n);
        info_set{i_user}(idxBhatta(1:k)) = 1;
    end
    
    for realization = 1:num_realizations
        if mod(realization, ceil(num_realizations/10)) == 1
            disp(['Sim iteration running = ', num2str(realization)]);
        end
        
        % encoding
        i_user = mod(realization, R) + 1;
        u = zeros(1, n);
        u(info_set{i_user} == 1) = randi([0 1], 1, k);
        x = polarEncode(u);
        
        % channel X_i -> Y_i is equivalent to BI-AWGN channel with input power equal to P_opt_user(i,i_p)
        y = sqrt(P_opt_total(i_user,i_p))*(1-2*double(x)) + randn(1,n);
        
        % decoding
        LLRin = 2*sqrt(P_opt_total(i_user,i_p))*y;
        
        if scl_flag == 0 || list_size == 1
            [uhat, ~, ~] = polarDecodeSSC(LLRin, info_set{i_user});
        else
            obj = PolarCode_SCL(n, k, info_set{i_user}, 0, []);
            [uhat, ~, ~] = obj.decode_scl_llr(LLRin', list_size);
            uhat = uhat';
        end
        
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
    Pe_bit(i_p) = nb_bit_err(i_p)/(realization*k);
    
    if (Pe_block(i_p) < 1e-4)
        break;
    end
end

toc
