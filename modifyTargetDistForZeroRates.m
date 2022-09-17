function output_dist = modifyTargetDistForZeroRates(input_dist, rate)
%function that modifies the target distribution when zero code rates are
%allocated

[num_dist, dist_size] = size(input_dist);
R = log2(dist_size);

if all(rate(:,2) == 0)
    output_dist = ones(size(input_dist))/dist_size;
elseif rate(1,2) == 0 && R == 2
    output_dist = zeros(num_dist, dist_size);
    for i_p = 1:num_dist
        beta = input_dist(i_p,2)/(input_dist(i_p,1) + input_dist(i_p,2));
        gamma = input_dist(i_p,3)/(input_dist(i_p,3) + input_dist(i_p,4));
        output_dist(i_p,:) = [0.5*(1-beta) 0.5*beta 0.5*gamma 0.5*(1-gamma)];
    end
elseif rate(2,2) == 0 && R == 2
    output_dist = zeros(num_dist, dist_size);
    for i_p = 1:num_dist
        alpha = input_dist(i_p,3) + input_dist(i_p,4);
        output_dist(i_p,:) = [0.5*(1-alpha) 0.5*(1-alpha) 0.5*alpha 0.5*alpha];
    end
else
    output_dist = input_dist;
end

end

