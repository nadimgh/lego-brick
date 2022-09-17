%testing script for a Marton coding scheme with optimal precoding over a two-user broadcast channel

clearvars;
close all;

n = 1024;
sum_rate = 1;
g = 0.9;
H = [1 g; g 1];
P_dB = -10:0.5:15;
scl_flag = 0;
list_size = 1;
num_realizations = 1e3;

P = 10.^(P_dB/10);
[R, T] = size(H);

rate_cell = getMaximumAchievableUserRates(H, P_dB, 'marton_opt');
rate_1 = getUserRates(sum_rate, rate_cell, 'marton_opt');

rate_2 = rate_1;
rate_2(1,2) = 0; rate_2(2,2) = rate_1(2,2) - 3*rate_1(2,2)/4; rate_2(2,1) = rate_1(2,1) - rate_1(1,2) - 3*rate_1(2,2)/4;

target_dist = rate_cell{1}.input_dist_marton_opt_precoding;
W_precoding = rate_cell{3};
x_all = 1-2*transpose(get_tuples(R));

nb_block_err = zeros(1, length(P));
nb_bit_err = zeros(1, length(P));
nb_block_err_user = zeros(R, length(P));
nb_bit_err_user = zeros(R, length(P));
Pe_block = zeros(1, length(P));
Pe_bit = zeros(1, length(P));
Pe_block_user = zeros(R, length(P));
Pe_bit_user = zeros(R, length(P));
true_dist = zeros(length(P), 2^R);

tic
for i_p = 1:length(P)
    if P_dB(i_p) < 7
        rate = rate_1;
    else
        rate = rate_2;
    end
    
    k = n*rate;
    target_dist(i_p,:) = modifyTargetDistForZeroRates(target_dist(i_p,:), rate);
    rho = target_dist(i_p,1) + target_dist(i_p,4) - (target_dist(i_p,2) + target_dist(i_p,3));
    rho_mat = [1 rho; rho 1];
    H_ch_effective = sqrt(P(i_p)/trace(rho_mat*(W_precoding{i_p}'*W_precoding{i_p}))) * H * W_precoding{i_p};
    
    gauss_params = cell(1,R);
    bhatta_pcm = cell(R,2);
    p_x_target = cell(1,R); % i-th element is p(x_1^i)
    p_x_user = cell(1,R); % i-th element is p(x_i)
    for i_user = 1:R
        idx0 = (x_all(i_user,:) == 1);
        mean0 = H_ch_effective(i_user,:) * x_all(:,idx0);
        alpha0 = target_dist(i_p,idx0)/sum(target_dist(i_p,idx0));
        sigma0 = ones(1, 2^(R-1));
        
        idx1 = (x_all(i_user,:) == -1);
        mean1 = H_ch_effective(i_user,:) * x_all(:,idx1);
        alpha1 = target_dist(i_p,idx1)/sum(target_dist(i_p,idx1));
        sigma1 = ones(1, 2^(R-1));
        
        p_x_user{i_user} = [sum(target_dist(i_p,idx0)) sum(target_dist(i_p,idx1))];
        
        gauss_params{i_user} = struct('p_x', p_x_user{i_user}, 'mean0', mean0, 'alpha0', alpha0, 'sigma0', sigma0, 'mean1', mean1, 'alpha1', alpha1, 'sigma1', sigma1);
        
        if i_user < R
            p_x_target{i_user} = sum(reshape(target_dist(i_p,:), 2^(R-i_user), []));
        else
            p_x_target{i_user} = target_dist(i_p,:);
        end
        
        bhatta_pcm{i_user,1} = computeBhattacharyyaPolarGaussian(n, p_x_user{i_user}, gauss_params{i_user}, 1);
        bhatta_pcm{i_user,2} = computeBhattacharyyaPolar(n, p_x_target{i_user});
    end
    
    H_struct = getParityCheckStruct(n, rate, bhatta_pcm);
    H_struct = StructToDouble(H_struct);
    
    for realization = 1:num_realizations
        if mod(realization, ceil(num_realizations/10)) == 1
            disp(['Sim iteration running = ', num2str(realization)]);
        end
        
        % Encoder
        u = cell(1,R);
        for i_user = 1:R
            u{i_user} = randi([0 1], k(i_user,1)-k(i_user,2), 1, 'int8');
        end
        [x, v] = MartonEncoder(u, target_dist(i_p,:), H_struct, 0, list_size);
        
        % Compute true input distribution
        x_mat = zeros(n, R);
        for i_user = 1:R
            x_mat(:,i_user) = x{i_user};
        end
        x_dec = bi2de(x_mat, 'left-msb')+1;
        [counts, indices] = groupcounts(x_dec);
        true_dist(i_p,indices) = true_dist(i_p,indices) + counts'/n;
        
        % Channel
        y = cell(1, R);
        for i_user = 1:R
            y{i_user} = zeros(n, 1);
            for i_antenna = 1:T
                y{i_user} = y{i_user} + H_ch_effective(i_user, i_antenna) * (1-2*double(x{i_antenna}));
            end
            y{i_user} = y{i_user} + randn(n,1);
        end
        
        % Decoder
        uhat = MartonDecoder(y, v, H_struct, gauss_params, scl_flag, list_size);
        
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
    
    Pe_block_user(:,i_p) = nb_block_err_user(:,i_p)/realization;
    Pe_block(i_p) = nb_block_err(i_p)/realization;
    Pe_bit_user(:,i_p) = nb_bit_err_user(:,i_p) ./ (realization * (k(:,1) - k(:,2)));
    Pe_bit(i_p) = nb_bit_err(i_p) / (realization * sum(k(:,1) - k(:,2)));
    true_dist(i_p,:) = true_dist(i_p,:)/realization;
    
    if (Pe_block(i_p) < 1e-4)
        break;
    end
end

toc

