function joint_dist = setQuantizationDistFromConditionalDistUplinkCRAN(input_dist, quantization_conditional_dist_user, params, power)
% set the joint input-quantization distribution of uplink C-RAN given the
% conditional and input distributions

e1 = zeros(4,1); e1(1) = 1;
e2 = zeros(4,1); e2(2) = 1;
e3 = zeros(4,1); e3(3) = 1;
e4 = zeros(4,1); e4(4) = 1;

input_dist_00 = input_dist(params)*e1;
input_dist_01 = input_dist(params)*e2;
input_dist_10 = input_dist(params)*e3;
input_dist_11 = input_dist(params)*e4;

quantization_conditional_dist_user1_00 = quantization_conditional_dist_user{1}(params,power)*e1;
quantization_conditional_dist_user1_01 = quantization_conditional_dist_user{1}(params,power)*e2;
quantization_conditional_dist_user1_10 = quantization_conditional_dist_user{1}(params,power)*e3;
quantization_conditional_dist_user1_11 = quantization_conditional_dist_user{1}(params,power)*e4;

quantization_conditional_dist_user2_00 = quantization_conditional_dist_user{2}(params,power)*e1;
quantization_conditional_dist_user2_01 = quantization_conditional_dist_user{2}(params,power)*e2;
quantization_conditional_dist_user2_10 = quantization_conditional_dist_user{2}(params,power)*e3;
quantization_conditional_dist_user2_11 = quantization_conditional_dist_user{2}(params,power)*e4;

joint_dist = [input_dist_00*(1-quantization_conditional_dist_user1_00)*(1-quantization_conditional_dist_user2_00) ...
    input_dist_00*(1-quantization_conditional_dist_user1_00)*quantization_conditional_dist_user2_00 ...
    input_dist_00*quantization_conditional_dist_user1_00*(1-quantization_conditional_dist_user2_00) ...
    input_dist_00*quantization_conditional_dist_user1_00*quantization_conditional_dist_user2_00 ...
    input_dist_01*(1-quantization_conditional_dist_user1_01)*(1-quantization_conditional_dist_user2_01) ...
    input_dist_01*(1-quantization_conditional_dist_user1_01)*quantization_conditional_dist_user2_01 ...
    input_dist_01*quantization_conditional_dist_user1_01*(1-quantization_conditional_dist_user2_01) ...
    input_dist_01*quantization_conditional_dist_user1_01*quantization_conditional_dist_user2_01 ...
    input_dist_10*(1-quantization_conditional_dist_user1_10)*(1-quantization_conditional_dist_user2_10) ...
    input_dist_10*(1-quantization_conditional_dist_user1_10)*quantization_conditional_dist_user2_10 ...
    input_dist_10*quantization_conditional_dist_user1_10*(1-quantization_conditional_dist_user2_10) ...
    input_dist_10*quantization_conditional_dist_user1_10*quantization_conditional_dist_user2_10 ...
    input_dist_11*(1-quantization_conditional_dist_user1_11)*(1-quantization_conditional_dist_user2_11) ...
    input_dist_11*(1-quantization_conditional_dist_user1_11)*quantization_conditional_dist_user2_11 ...
    input_dist_11*quantization_conditional_dist_user1_11*(1-quantization_conditional_dist_user2_11) ...
    input_dist_11*quantization_conditional_dist_user1_11*quantization_conditional_dist_user2_11];
joint_dist = joint_dist/sum(joint_dist);

end

