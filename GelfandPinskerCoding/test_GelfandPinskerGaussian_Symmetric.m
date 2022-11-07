%testing script for the symmetric coding scheme over the Gaussian channel Y=X+gS+Z, where X and S are binary and belong to {+/-sqrt(P)}

clearvars;
close all;

current_folder = pwd;
addpath(genpath([fileparts(current_folder), filesep, 'src' ]));

n = 1024;
rate = 1/2;
P_dB = -10:0.5:15;
g = 0.9;
theta = 0.1;
scl_flag = 0;
list_size = 1;
num_realizations = 5e4;

p_s = [1-theta theta];
k = round(n*rate); k1 = k(1);
P = 10.^(P_dB/10);
theta = p_s(2);

design_dist = [0.5*(1-theta) 0.5*theta 0.5*theta 0.5*(1-theta)];
p_x = [0.5 0.5];
sigma00 = 1; sigma01 = 1; sigma10 = 1; sigma11 = 1;
alpha00 = 1; alpha01 = 1; alpha10 = 1; alpha11 = 1;

nb_block_err = zeros(1, length(P));
nb_bit_err = zeros(1, length(P));
Pe_block = zeros(1, length(P));
Pe_bit = zeros(1, length(P));

tic
for i_p = 1:length(P)
    mean00 = sqrt(P(i_p))*(1+g); mean01 = sqrt(P(i_p))*(-1+g); mean10 = sqrt(P(i_p))*(1-g); mean11 = sqrt(P(i_p))*(-1-g);
    gauss_params = struct('p_x', design_dist, 'mean00', mean00, 'mean01', mean01, 'mean10', mean10, 'mean11', mean11, 'sigma00', sigma00, 'sigma01', sigma01, 'sigma10', sigma10, 'sigma11', sigma11, ...
        'alpha00', alpha00, 'alpha01', alpha01, 'alpha10', alpha10, 'alpha11', alpha11);
    gauss_params_XY = struct('p_x', p_x, 'mean0', [mean00 mean10], 'sigma0', [sigma00 sigma10], 'alpha0', [(1-theta)*alpha00 theta*alpha10], ...
        'mean1', [mean01 mean11], 'sigma1', [sigma01 sigma11], 'alpha1', [(1-theta)*alpha01 theta*alpha11]);
    
    bhatta = computeBhattacharyyaPolarGaussian(n, p_x, gauss_params_XY, 0);
    [~, idxBhatta] = sort(bhatta);
    info_set = zeros(1,n);
    info_set(idxBhatta(1:k1)) = 1;

    for realization = 1:num_realizations
        if mod(realization, ceil(num_realizations/10)) == 1
            disp(['Sim iteration running = ', num2str(realization)]);
        end

        s = randsmpl(p_s, 1, n, 'int8')-1;    % independent sampling from p_s
        
        u = zeros(1, n);
        u(info_set == 1) = randi([0 1], 1, k1);

        x = polarEncode(u);

        y = generateGaussianMixtureChannelOutput(x, s, gauss_params);
        
        summ0 = zeros(1,n);
        summ1 = zeros(1,n);
        for i = 1:length(gauss_params_XY.mean0)
            summ0 = summ0 + gauss_params_XY.alpha0(i)*normpdf(y, gauss_params_XY.mean0(i), gauss_params_XY.sigma0(i));
            summ1 = summ1 + gauss_params_XY.alpha1(i)*normpdf(y, gauss_params_XY.mean1(i), gauss_params_XY.sigma1(i));
        end
        LLRin = log(summ0./summ1);
        
        if scl_flag == 0 || list_size == 1
            [uhat, ~, ~] = polarDecodeSSC(LLRin, info_set);
        else
            obj = PolarCode_SCL(n, k1, info_set, 0, []);
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
    Pe_bit(i_p) = nb_bit_err(i_p)/(realization*k1);
    
    if (Pe_block(i_p) < 1e-4)
        break;
    end
end

toc

