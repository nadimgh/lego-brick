%compute the maximum achievable rates for the downlink C-RAN problem using the genetic algorithm

g = 0.9;
H = [1 g; g 1];
P_dB = -10;%:0.5:15;
C_vec = [1 1];
[R, T] = size(H);
P = 10.^(P_dB/10);

rng('shuffle');
tol = 1e-5;

s_all = 1-2*transpose(get_tuples(R+T));
input_dist = @(params) setInputDistFromParamsCRAN(params);
W_mat = @(w) [sqrt(w(1)) 0; 0 sqrt(w(2))];
x_all = @(w) W_mat(w) * s_all(R+1:R+T,:);
low_val = @(w) min(H * W_mat(w) * s_all(R+1:R+T,:), [], 'all') - 5;
upper_val = @(w) max(H * W_mat(w) * s_all(R+1:R+T,:), [], 'all') + 5;
sum_rate_fun = @(params,w) 0;
for i_user = 1:R
    mean_vec = @(w) H(i_user,:) * x_all(w);
    sigma_vec = ones(1,2^(R+T));
    alpha_vec = input_dist;
    fun = entropy_function(mean_vec, sigma_vec, alpha_vec);
    entropy_Y = @(params,w) integral(@(y) fun(y,params,w), low_val(w), upper_val(w));

    idx_plus = (s_all(i_user,:) == 1);
    matrix_plus = zeros(2^(R+T), 2^(R+T-1)); matrix_plus(idx_plus,:) = eye(2^(R+T-1));
    mean_vec_plus = @(w) H(i_user,:) * x_all(w) * matrix_plus;
    sigma_vec_plus = ones(1,2^(R+T-1));
    alpha_vec_plus = @(params) (input_dist(params)*matrix_plus)/sum(input_dist(params)*matrix_plus);
    low_val_plus = @(w) min(mean_vec_plus(w), [], 'all') - 5;
    upper_val_plus = @(w) max(mean_vec_plus(w), [], 'all') + 5;
    fun_plus = entropy_function(mean_vec_plus, sigma_vec_plus, alpha_vec_plus);
    idx_minus = (s_all(i_user,:) == -1);
    matrix_minus = zeros(2^(R+T), 2^(R+T-1)); matrix_minus(idx_minus,:) = eye(2^(R+T-1));
    mean_vec_minus = @(w) H(i_user,:) * x_all(w) * matrix_minus;
    sigma_vec_minus = ones(1,2^(R+T-1));
    alpha_vec_minus = @(params) (input_dist(params)*matrix_minus)/sum(input_dist(params)*matrix_minus);
    low_val_minus = @(w) min(mean_vec_minus(w), [], 'all') - 5;
    upper_val_minus = @(w) max(mean_vec_minus(w), [], 'all') + 5;
    fun_minus = entropy_function(mean_vec_minus, sigma_vec_minus, alpha_vec_minus);
    entropy_YgivenX = @(params,w) sum(input_dist(params)*matrix_plus) * integral(@(y) fun_plus(y,params,w), low_val_plus(w), upper_val_plus(w)) + ...
        sum(input_dist(params)*matrix_minus) * integral(@(y) fun_minus(y,params,w), low_val_minus(w), upper_val_minus(w));

    if i_user == 1
        sum_rate_fun = @(params,w) sum_rate_fun(params,w) + entropy_Y(params,w) - entropy_YgivenX(params,w);
    else
        input_dist_user = @(params) sum(reshape(input_dist(params), 2^(R+T-i_user), []));
        sum_rate_fun = @(params,w) sum_rate_fun(params,w) + entropy_Y(params,w) - entropy_YgivenX(params,w) - computeMutualInformation(input_dist_user(params));
    end
end

tic
params = zeros(length(P), 2^(R+T)-1);
input_dist = zeros(length(P), 2^(R+T));
sum_rate = zeros(1, length(P));
W_opt = cell(1,length(P));
for i_p = 1:length(P)
    disp(['Sim iteration running = ', num2str(i_p)]);
    
    options = optimoptions('ga', 'PopulationSize', 100, 'HybridFcn', @fmincon, 'ConstraintTolerance', tol);
    
    fun_obj = @(x) -1*sum_rate_fun(x(1:2^(R+T)-1), x(2^(R+T):end));
    [x_opt, fval] = ga(fun_obj, 2^(R+T)+R-1, [], [], [zeros(1,2^(R+T)-1) 1 1], P(i_p), tol*ones(1,2^(R+T)+R-1), [ones(1,2^(R+T)-1)-tol P(i_p) P(i_p)], @(x)computeNonLinearConstraintFunction(x, C_vec), options);
    
    sum_rate(i_p) = -fval;
    params(i_p,:) = x_opt(1:2^(R+T)-1);
    input_dist(i_p,:) = setInputDistFromParamsCRAN(params(i_p,:));
    W_opt{i_p} = [sqrt(x_opt(2^(R+T))) 0; 0 sqrt(x_opt(2^(R+T)+1))];
end

matfile = sprintf('AchievableRates_DownlinkCRAN_GA_g%.2f_C%.2f.mat', g, C_vec(1));
if ~exist(matfile, 'file')
    save(matfile, 'params', 'input_dist', 'sum_rate', 'W_opt', 'P_dB', 'H', 'C_vec');
end
toc



function z = pdf_channel_output(y, params, w, mean_vec, sigma_vec, alpha_vec)
    z = 0 ;
    mean_temp = mean_vec(w);
    alpha_temp = alpha_vec(params);
    for k = 1:length(sigma_vec)
        z = z + alpha_temp(k) .* 1/sqrt(2*pi*sigma_vec(k)^2) .* exp( -(y-mean_temp(k)).^2/(2*sigma_vec(k)^2) );
    end
end

function z = entropy_function(mean_vec, sigma_vec, alpha_vec)
    z = @(y,params,w) max(0, -1 * pdf_channel_output(y, params, w, mean_vec, sigma_vec, alpha_vec) .* log2(pdf_channel_output(y, params, w, mean_vec, sigma_vec, alpha_vec)) );
end

