%compute the maximum achievable sum-rates of different coding stratgies
%for the Gaussian broadcast channel with channel gain matrix H and sum-power constraints P_dB.
%Coding strategies considered: 
%  1- Marton coding with optimal precoding,
%  2- Marton coding (no precoding), 
%  3- Symmetric coding with optimal precoding,
%  4- Symmetric coding (no precoding), 
%  5- MMSE precoding, 
%  6- Zero-forcing precoding, 
%  7- Time-division

g = 0.9;
H = [1 g; g 1];
P_dB = -10:0.5:15;
[R, T] = size(H);

P = 10.^(P_dB/10);

input_dist_marton = zeros(length(P), 4);
input_dist_marton_opt_precoding = zeros(length(P), 4);
rate_marton = zeros(1, length(P));
rate_marton_opt_precoding = zeros(1, length(P));
rate_symmetric = zeros(1, length(P));
rate_symmetric_opt_precoding = zeros(1, length(P));
rate_mmse = zeros(1, length(P));
rate_time_division = zeros(1, length(P));
rate_zf = zeros(1, length(P));
W_opt = cell(1, length(P));
W_opt_symmetric = cell(1, length(P));

tic
for i_p = 1:length(P)
    disp(['Sim iteration running = ', num2str(i_p)]);
    
    [rate_marton_opt_precoding(i_p), input_dist_marton_opt_precoding(i_p,:), W_opt{i_p}] = optimizeInputDistAndPrecodingMatrix_PSO(H, P(i_p), 'marton_opt');
    
    [rate_marton(i_p), input_dist_marton(i_p,:), ~] = optimizeInputDistAndPrecodingMatrix_PSO(H, P(i_p), 'marton');
    
    [rate_symmetric_opt_precoding(i_p), ~, W_opt_symmetric{i_p}] = optimizeInputDistAndPrecodingMatrix_PSO(H, P(i_p), 'symmetric_opt');
    
    [rate_symmetric(i_p), ~] = computeCapacityMartonSpecificInputDist(H, P(i_p), eye(R), ones(1,2^R)/(2^R));
    
    W_mmse = (H'*H + (R/P(i_p))*eye(T))\(H');
    [rate_mmse(i_p), ~] = computeCapacityPrecoding(H, P(i_p), W_mmse);
    
    [rate_zf(i_p), ~] = computeCapacityZeroForcing(H, P(i_p));
    
    rate_time_division(i_p) = computeCapacityTimeDivision(H, P(i_p));
end

matfile = sprintf('AchievableRates_BC_g%.2f.mat', g);
if ~exist(matfile, 'file')
    save(matfile, 'rate_marton', 'rate_marton_opt_precoding', 'rate_symmetric', 'rate_symmetric_opt_precoding', ...
        'input_dist_marton_opt_precoding', 'input_dist_marton', 'W_opt', 'W_opt_symmetric', 'rate_mmse', 'rate_zf', 'rate_time_division', 'g', 'P_dB', 'H');
end
toc
