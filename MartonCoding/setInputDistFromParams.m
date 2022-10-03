function input_dist = setInputDistFromParams(params)
%get the input distribution from the parameters of the optimization problem
%corresponding to the capacity computation of the CRAN problem

input_dist = [(1-params(1))*(1-params(2)) (1-params(1))*params(2) params(1)*(1-params(3)) params(1)*params(3)];
input_dist = input_dist/sum(input_dist);

end

