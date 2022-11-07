function [c, ceq] = computeNonLinearConstraintFunction(x, C_vec)
%function of non-linear constraints to be used in optimization problem
%(corresponds to the fact the input distribution should satisfy compression constraints of C-RAN problem)

params = x(1:15);
input_dist = setInputDistFromParamsDownlinkCRAN(params);
input_dist_u1u2x1 = sum(reshape(input_dist, 2, []));
constraint_1 = computeMutualInformation(input_dist_u1u2x1);
constraint_2 = computeMutualInformation(input_dist);

c(1) = constraint_1 - C_vec(1);
c(2) = constraint_2 - C_vec(2);
ceq = [];

end
