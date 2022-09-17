%testing script for MMSE precoding scheme over a two-user broadcast channel

clearvars;
close all;

n = 1024;
sum_rate = 1;
g = 0.9;
H = [1 g; g 1];
P_dB = -10:0.5:25;
scl_flag = 0;
list_size = 1;
num_realizations = 1e3;

P = 10.^(P_dB/10);
[R, T] = size(H);

rate_cell = getMaximumAchievableUserRates(H, P_dB, 'mmse');
rate = getUserRates(sum_rate, rate_cell, 'mmse');

W_precoding = cell(1, length(P));
for i_p = 1:length(P)
    W_precoding{i_p} = (H'*H + (R/P(i_p))*eye(T))\(H'); % MMSE precoding matrix
end

x_all = 1-2*transpose(get_tuples(R));

nb_block_err = zeros(1, length(P));
nb_bit_err = zeros(1, length(P));
nb_block_err_user = zeros(R, length(P));
nb_bit_err_user = zeros(R, length(P));
Pe_block = zeros(1, length(P));
Pe_bit = zeros(1, length(P));
Pe_block_user = zeros(R, length(P));
Pe_bit_user = zeros(R, length(P));

tic
for i_p = 1:length(P)
    H_ch_effective = H*W_precoding{i_p};
    k = n*rate;
    
    gauss_params = cell(1,R);
    info_set = cell(1,R);
    for i_user = 1:R
        idx0 = (x_all(i_user,:) == 1);
        mean0 = (sqrt(P(i_p))/norm(W_precoding{i_p}(:))) * H_ch_effective(i_user,:) * x_all(:,idx0);
        alpha0 = ones(1, 2^(R-1))/(2^(R-1));
        sigma0 = ones(1, 2^(R-1));
        
        idx1 = (x_all(i_user,:) == -1);
        mean1 = (sqrt(P(i_p))/norm(W_precoding{i_p}(:))) * H_ch_effective(i_user,:) * x_all(:,idx1);
        alpha1 = ones(1, 2^(R-1))/(2^(R-1));
        sigma1 = ones(1, 2^(R-1));
        
        gauss_params{i_user} = struct('p_x', [0.5 0.5], 'mean0', mean0, 'alpha0', alpha0, 'sigma0', sigma0, 'mean1', mean1, 'alpha1', alpha1, 'sigma1', sigma1);
        
        bhatta = computeBhattacharyyaPolarGaussian(n, [0.5 0.5], gauss_params{i_user}, 0);
        [~, idxBhatta] = sort(bhatta);
        info_set{i_user} = zeros(1,n);
        info_set{i_user}(idxBhatta(1:k(i_user))) = 1;
    end
    
    for realization = 1:num_realizations
        if mod(realization, ceil(num_realizations/10)) == 1
            disp(['Sim iteration running = ', num2str(realization)]);
        end
        
        % initialize
        u = cell(1, R);
        x = cell(1, R);
        uhat = cell(1, R);
        err = zeros(1, R);
        
        % encoding
        for i_user = 1:R
            u{i_user} = zeros(1, n);
            u{i_user}(info_set{i_user} == 1) = randi([0 1], 1, k(i_user));
            x{i_user} = polarEncode(u{i_user});
        end
        
        % Channel
        y = cell(1, R);
        for i_user = 1:R
            y{i_user} = zeros(n, 1);
            for i_antenna = 1:T
                y{i_user} = y{i_user} + (sqrt(P(i_p))/norm(W_precoding{i_p}(:))) * H_ch_effective(i_user, i_antenna) * (1-2*double(x{i_antenna}'));
            end
            y{i_user} = y{i_user} + randn(n,1);
            y{i_user} = y{i_user}';
        end
        
        % decoding
        for i_user = 1:R
            summ0 = zeros(1,n);
            summ1 = zeros(1,n);
            for i = 1:length(gauss_params{i_user}.mean0)
                summ0 = summ0 + gauss_params{i_user}.alpha0(i)*normpdf(y{i_user}, gauss_params{i_user}.mean0(i), gauss_params{i_user}.sigma0(i));
                summ1 = summ1 + gauss_params{i_user}.alpha1(i)*normpdf(y{i_user}, gauss_params{i_user}.mean1(i), gauss_params{i_user}.sigma1(i));
            end
            LLRin = log(summ0./summ1);
            
            if scl_flag == 0 || list_size == 1
                [uhat{i_user}, ~, ~] = polarDecodeSSC(LLRin, info_set{i_user});
            else
                obj = PolarCode_SCL(n, k(i_user), info_set{i_user}, 0, []);
                [uhat{i_user}, ~, ~] = obj.decode_scl_llr(LLRin', list_size);
                uhat{i_user} = uhat{i_user}';
            end
            
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
    Pe_bit_user(:,i_p) = nb_bit_err_user(:,i_p) ./ (realization * k);
    Pe_bit(i_p) = nb_bit_err(i_p) / (realization * sum(k));
    
    if (Pe_block(i_p) < 1e-4)
        break;
    end
end

toc
