function input_dist = setInputDistFromParamsUplinkCRAN(params)
%compute the input distribution from the parameters of the optimization problem
%corresponding to the capacity computation of the uplink CRAN problem

input_dist = [(1-params(1))*(1-params(2)) (1-params(1))*params(2) params(1)*(1-params(2)) params(1)*params(2)];
input_dist = input_dist/sum(input_dist);

end

