function [sum_rate, input_dist, W_opt] = optimizeInputDistAndPrecodingMatrix_PSO(H, P, scheme)
%compute the maximum achievable rates for Marton coding and symmetric coding using particle swarm optimization

[R, ~] = size(H);

rng('shuffle');

s_all = 1-2*transpose(get_tuples(R));

input_dist = @(params) setInputDistFromParams(params);
rho = @(params) (1-params(1))*(1-params(2)) + params(1)*params(3) - ( (1-params(1))*params(2) + params(1)*(1-params(3)) );
rho_mat = @(params) [1 rho(params); rho(params) 1];
W_mat = @(w) [w(1) w(3); w(2) w(4)];
low_val = @(power,params,w) min(sqrt( power/trace(rho_mat(params)*W_mat(w)'*W_mat(w)) ) * H * W_mat(w) * s_all, [], 'all') - 5;
upper_val = @(power,params,w) max(sqrt( power/trace(rho_mat(params)*W_mat(w)'*W_mat(w)) ) * H * W_mat(w) * s_all, [], 'all') + 5;
sum_rate_fun = @(power,params,w) 0;
for i_user = 1:R
    mean_vec = @(power,params,w) sqrt( power/trace(rho_mat(params)*W_mat(w)'*W_mat(w)) ) * H(i_user,:) * W_mat(w) * s_all;
    sigma_vec = ones(1,2^R);
    alpha_vec = input_dist;
    fun = entropy_function(mean_vec, sigma_vec, alpha_vec);
    entropy_Y = @(power,params,w) integral(@(y) fun(y,power,params,w), low_val(power,params,w), upper_val(power,params,w));

    idx_plus = (s_all(i_user,:) == 1);
    mean_vec_plus = @(power,params,w) sqrt( power/trace(rho_mat(params)*W_mat(w)'*W_mat(w)) ) * H(i_user,:) * W_mat(w) * s_all(:,idx_plus);
    sigma_vec_plus = ones(1,2^(R-1));
    matrix_plus = zeros(2^R, 2^(R-1)); matrix_plus(idx_plus,:) = eye(2^(R-1));
    alpha_vec_plus = @(params) (input_dist(params)*matrix_plus)/sum(input_dist(params)*matrix_plus);
    low_val_plus = @(power,params,w) min(mean_vec_plus(power,params,w), [], 'all') - 5;
    upper_val_plus = @(power,params,w) max(mean_vec_plus(power,params,w), [], 'all') + 5;
    fun_plus = entropy_function(mean_vec_plus, sigma_vec_plus, alpha_vec_plus);
    idx_minus = (s_all(i_user,:) == -1);
    mean_vec_minus = @(power,params,w) sqrt( power/trace(rho_mat(params)*W_mat(w)'*W_mat(w)) ) * H(i_user,:) * W_mat(w) * s_all(:,idx_minus);
    sigma_vec_minus = ones(1,2^(R-1));
    matrix_minus = zeros(2^R, 2^(R-1)); matrix_minus(idx_minus,:) = eye(2^(R-1));
    alpha_vec_minus = @(params) (input_dist(params)*matrix_minus)/sum(input_dist(params)*matrix_minus);
    low_val_minus = @(power,params,w) min(mean_vec_minus(power,params,w), [], 'all') - 5;
    upper_val_minus = @(power,params,w) max(mean_vec_minus(power,params,w), [], 'all') + 5;
    fun_minus = entropy_function(mean_vec_minus, sigma_vec_minus, alpha_vec_minus);
    entropy_YgivenX = @(power,params,w) sum(input_dist(params)*matrix_plus) * integral(@(y) fun_plus(y,power,params,w), low_val_plus(power,params,w), upper_val_plus(power,params,w)) + ...
        sum(input_dist(params)*matrix_minus) * integral(@(y) fun_minus(y,power,params,w), low_val_minus(power,params,w), upper_val_minus(power,params,w));

    if i_user == 1
        sum_rate_fun = @(power,params,w) sum_rate_fun(power,params,w) + entropy_Y(power,params,w) - entropy_YgivenX(power,params,w);
    elseif i_user < R
        input_dist_user = @(params) sum(reshape(input_dist(params), 2^(R-i_user), []));
        sum_rate_fun = @(power,params,w) sum_rate_fun(power,params,w) + entropy_Y(power,params,w) - entropy_YgivenX(power,params,w) - computeMutualInformation(input_dist_user(params));
    else
        sum_rate_fun = @(power,params,w) sum_rate_fun(power,params,w) + entropy_Y(power,params,w) - entropy_YgivenX(power,params,w) - computeMutualInformation(input_dist(params));
    end
end

params = zeros(length(P), 2^R-1);
input_dist = zeros(length(P), 2^R);
sum_rate = zeros(1, length(P));
W_opt = cell(1,length(P));
max_val = 2;

for i_p = 1:length(P)    
    if strcmp(scheme, 'marton_opt')
        if i_p > 1
            options = optimoptions('particleswarm', 'SwarmSize', 200, 'HybridFcn', @fmincon, 'InitialSwarmMatrix', [params(i_p-1,:) transpose(W_opt{i_p-1}(:))]);
        else
            options = optimoptions('particleswarm', 'SwarmSize', 200, 'HybridFcn', @fmincon);
        end
        fun_obj = @(x) -1*sum_rate_fun(P(i_p), x(1:3), x(4:7));
        [x_opt, fval] = particleswarm(fun_obj, 7, [zeros(1,3) -max_val*ones(1,4)], [ones(1,3) max_val*ones(1,4)], options);
        
        sum_rate(i_p) = -fval;
        params(i_p,:) = x_opt(1:3);
        input_dist(i_p,:) = setInputDistFromParams(params(i_p,:));
        rho = input_dist(i_p,1) + input_dist(i_p,4) - (input_dist(i_p,2) + input_dist(i_p,3));
        W_opt{i_p} = reshape(x_opt(4:7), R, []);
        W_opt{i_p} = W_opt{i_p}/sqrt(trace([1 rho; rho 1]*W_opt{i_p}'*W_opt{i_p}));
    elseif strcmp(scheme, 'marton')
        if i_p > 1
            options = optimoptions('particleswarm', 'SwarmSize', 200, 'HybridFcn', @fmincon, 'InitialSwarmMatrix', params(i_p-1,:));
        else
            options = optimoptions('particleswarm', 'SwarmSize', 200, 'HybridFcn', @fmincon);
        end
        fun_obj = @(x) -1*sum_rate_fun(P(i_p), x(1:3), [1 0 0 1]/sqrt(2));
        [x_opt, fval] = particleswarm(fun_obj, 3, zeros(1,3), ones(1,3), options);
        
        sum_rate(i_p) = -fval;
        params(i_p,:) = x_opt(1:3);
        input_dist(i_p,:) = setInputDistFromParams(params(i_p,:));
        W_opt{i_p} = [1 0; 0 1]/sqrt(2);
    elseif strcmp(scheme, 'symmetric_opt')
        if i_p > 1
            options = optimoptions('particleswarm', 'SwarmSize', 200, 'HybridFcn', @fmincon, 'InitialSwarmMatrix', transpose(W_opt{i_p-1}(:)));
        else
            options = optimoptions('particleswarm', 'SwarmSize', 200, 'HybridFcn', @fmincon);
        end
        fun_obj = @(x) -1*sum_rate_fun(P(i_p), [0.5 0.5 0.5], x(1:4));
        [x_opt, fval] = particleswarm(fun_obj, 4, -max_val*ones(1,4), max_val*ones(1,4), options);
        
        sum_rate(i_p) = -fval;
        params(i_p,:) = [0.5 0.5 0.5];
        input_dist(i_p,:) = setInputDistFromParams(params(i_p,:));
        W_opt{i_p} = reshape(x_opt(1:4), R, []);
        W_opt{i_p} = W_opt{i_p}/sqrt(trace([1 0; 0 1]*W_opt{i_p}'*W_opt{i_p}));
    end
end

end




function z = pdf_channel_output(y, power, params, w, mean_vec, sigma_vec, alpha_vec)
    z = 0 ;
    mean_temp = mean_vec(power, params, w);
    alpha_temp = alpha_vec(params);
    for k = 1:length(sigma_vec)
        z = z + alpha_temp(k) .* 1/sqrt(2*pi*sigma_vec(k)^2) .* exp( -(y-mean_temp(k)).^2/(2*sigma_vec(k)^2) );
    end
end

function z = entropy_function(mean_vec, sigma_vec, alpha_vec)
    z = @(y,power,params,w) max(0, -1 * pdf_channel_output(y, power, params, w, mean_vec, sigma_vec, alpha_vec) .* log2(pdf_channel_output(y, power, params, w, mean_vec, sigma_vec, alpha_vec)) );
end

