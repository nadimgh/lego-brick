function [c, ceq] = computeNonLinearConstraintFunction(x, C_vec)
%function of non-linear constraints to be used in the genetic algorithm

params = x(1:15);
input_dist = setInputDistFromParamsCRAN(params);
input_dist_u1u2x1 = sum(reshape(input_dist, 2, []));
constraint_1 = computeMutualInformation(input_dist_u1u2x1);
constraint_2 = computeMutualInformation(input_dist);

c(1) = constraint_1 - C_vec(1); % replace with the value of the fronthaul constraints at each run
c(2) = constraint_2 - C_vec(2);
ceq = [];

end
