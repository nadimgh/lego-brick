function mutual_info = computeMutualInformation(joint_dist)
%compute the mutual information I(Y1,Y2;X) for the joint distribution
%joint_dist which takes the form p(y1,y2,x) (in this order)

if ( any(joint_dist < 0) || any(joint_dist > 1) || abs(sum(joint_dist) - 1) > 1e-5 )
    error('Joint_dist should be a probability distribution.');
end
joint_dist = joint_dist(:);

m = log2(length(joint_dist));
if (m ~= floor(m))
    error('Length of joint_dist should be a power of 2.');
end

p_x = sum(reshape(joint_dist, 2, []), 2);
p_y = transpose(sum(reshape(joint_dist, 2, [])));
product_dist = kron(p_y, p_x);

idx = (joint_dist ~= 0) & (product_dist ~= 0);
mutual_info = sum(joint_dist(idx) .* log2(joint_dist(idx)./product_dist(idx)));

end

