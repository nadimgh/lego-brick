function [c, ceq] = computeNonLinearConstraintFunction(x, H, P, C_vec)
%function of non-linear constraints to be in optimization solver
%(corresponds to the fact the target distributions should satisfy compression constraints of C-RAN problem)

[R, ~] = size(H);
gauss_params = cell(1,R);
quantization_conditional_dist_user = cell(1,R);
output_conditional_dist_user = cell(1,R);
s_all = 1-2*transpose(get_tuples(R));
input_dist = setInputDistFromParamsUplinkCRAN(x);
W_mat = [sqrt(P) 0; 0 sqrt(P)];
x_all = W_mat * s_all;
for i_user = 1:R
    mean_vec = H(i_user,:) * x_all;
    sigma_vec = ones(1,2^R);
    alpha_vec = input_dist;
    gauss_params{i_user} = struct('alpha_vec', alpha_vec, 'mean_vec', mean_vec, 'sigma_vec', sigma_vec);
    b = getLloydMaxQuantizationParameters(gauss_params{i_user});
    output_conditional_dist_user{i_user} = 1 - Q((b-mean_vec)./sigma_vec);
    quantization_conditional_dist_user{i_user} = x(2*i_user+1)*(1-output_conditional_dist_user{i_user}) + (1-x(2*i_user+2))*output_conditional_dist_user{i_user};
end

output_dist_eval = setOutputDistFromConditionalDistUplinkCRAN(input_dist, output_conditional_dist_user, x, P);
joint_output_quantization_dist_eval = setJointOutputQuantizationDistUplinkCRAN(output_dist_eval, x, P);
quantization_dist_eval = transpose(sum(reshape(joint_output_quantization_dist_eval, 4, []), 2));

dist_y1u1 = [sum(joint_output_quantization_dist_eval([1 2 5 6])) sum(joint_output_quantization_dist_eval([3 4 7 8])) ...
    sum(joint_output_quantization_dist_eval([9 10 13 14])) sum(joint_output_quantization_dist_eval([11 12 15 16]))];
dist_y2u2 = [sum(joint_output_quantization_dist_eval([1 3 9 11])) sum(joint_output_quantization_dist_eval([2 4 10 12])) ...
        sum(joint_output_quantization_dist_eval([5 7 13 15])) sum(joint_output_quantization_dist_eval([6 8 14 16]))];

constraint_1 = computeMutualInformation(dist_y1u1);
constraint_2 = computeMutualInformation(dist_y2u2) - computeMutualInformation(quantization_dist_eval);

c(1) = constraint_1 - C_vec(1); % replace with the value of the backhaul constraints at each run
c(2) = constraint_2 - C_vec(2);
ceq = [];

end

