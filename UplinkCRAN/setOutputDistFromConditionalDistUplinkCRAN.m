function output_dist = setOutputDistFromConditionalDistUplinkCRAN(input_dist, output_conditional_dist_user, params, power)
% set the channel output distribution of uplink C-RAN given the conditional and input distributions

if isa(input_dist, 'function_handle')
    input_dist = input_dist(params);
end
if isa(output_conditional_dist_user{1}, 'function_handle')
    output_conditional_dist_user{1} = output_conditional_dist_user{1}(params,power);
end
if isa(output_conditional_dist_user{2}, 'function_handle')
    output_conditional_dist_user{2} = output_conditional_dist_user{2}(params,power);
end

output_dist_00 = input_dist * ((1-output_conditional_dist_user{1}') .* (1-output_conditional_dist_user{2}'));
output_dist_01 = input_dist * ((1-output_conditional_dist_user{1}') .* (output_conditional_dist_user{2}'));
output_dist_10 = input_dist * ((output_conditional_dist_user{1}') .* (1-output_conditional_dist_user{2}'));
output_dist_11 = input_dist * ((output_conditional_dist_user{1}') .* (output_conditional_dist_user{2}'));

output_dist = [output_dist_00 output_dist_01 output_dist_10 output_dist_11];

end

