function joint_dist = setJointOutputQuantizationDistUplinkCRAN(output_dist, params, power)
% set the joint output-quantization distribution of uplink C-RAN given the
% channel output distribution

if isa(output_dist, 'function_handle')
    output_dist = output_dist(params,power);
end

e1 = zeros(4,1); e1(1) = 1;
e2 = zeros(4,1); e2(2) = 1;
e3 = zeros(4,1); e3(3) = 1;
e4 = zeros(4,1); e4(4) = 1;

output_dist_00 = output_dist*e1;
output_dist_01 = output_dist*e2;
output_dist_10 = output_dist*e3;
output_dist_11 = output_dist*e4;

% joint_dist = [P(Y1=0,Y2=0,Y1hat=0,Y2hat=0) P(Y1=0,Y2=0,Y1hat=0,Y2hat=1) P(Y1=0,Y2=0,Y1hat=1,Y2hat=0) P(Y1=0,Y2=0,Y1hat=1,Y2hat=1) ...
%               P(Y1=0,Y2=1,Y1hat=0,Y2hat=0) P(Y1=0,Y2=1,Y1hat=0,Y2hat=1) P(Y1=0,Y2=1,Y1hat=1,Y2hat=0) P(Y1=0,Y2=1,Y1hat=1,Y2hat=1) ...
%               P(Y1=1,Y2=0,Y1hat=0,Y2hat=0) P(Y1=1,Y2=0,Y1hat=0,Y2hat=1) P(Y1=1,Y2=0,Y1hat=1,Y2hat=0) P(Y1=1,Y2=0,Y1hat=1,Y2hat=1) ...
%               P(Y1=1,Y2=1,Y1hat=0,Y2hat=0) P(Y1=1,Y2=1,Y1hat=0,Y2hat=1) P(Y1=1,Y2=1,Y1hat=1,Y2hat=0) P(Y1=1,Y2=1,Y1hat=1,Y2hat=1) ]
joint_dist = [output_dist_00*(1-params(3))*(1-params(5)) ...
    output_dist_00*(1-params(3))*params(5) ...
    output_dist_00*params(3)*(1-params(5)) ...
    output_dist_00*params(3)*params(5) ...
    output_dist_01*(1-params(3))*params(6) ...
    output_dist_01*(1-params(3))*(1-params(6)) ...
    output_dist_01*params(3)*params(6) ...
    output_dist_01*params(3)*(1-params(6)) ...
    output_dist_10*params(4)*(1-params(5)) ...
    output_dist_10*params(4)*params(5) ...
    output_dist_10*(1-params(4))*(1-params(5)) ...
    output_dist_10*(1-params(4))*params(5) ...
    output_dist_11*params(4)*params(6) ...
    output_dist_11*params(4)*(1-params(6)) ...
    output_dist_11*(1-params(4))*params(6) ...
    output_dist_11*(1-params(4))*(1-params(6))];
joint_dist = joint_dist/sum(joint_dist);

end

