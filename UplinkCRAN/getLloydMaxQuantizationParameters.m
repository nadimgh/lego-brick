function [partition, codebook, output_dist] = getLloydMaxQuantizationParameters(gauss_params)
%find quantization symbols and boundaries using the Lloyd-Max algorithm

tol = 1e-5;
max_iterations = 1e3;

mean_vec = gauss_params.mean_vec;
sigma_vec = gauss_params.sigma_vec;
alpha_vec = gauss_params.alpha_vec;

b = 0;
mse_old = realmax;
mse_new = 100;
iter = 1;

while abs(mse_old - mse_new) > tol && iter <= max_iterations
    num = integral(mean_func(mean_vec, sigma_vec, alpha_vec), -Inf, b);
    den = integral(pdf_func(mean_vec, sigma_vec, alpha_vec), -Inf, b);
    a1 = num/den;
    
    num = integral(mean_func(mean_vec, sigma_vec, alpha_vec), b, Inf);
    den = integral(pdf_func(mean_vec, sigma_vec, alpha_vec), b, Inf);
    a2 = num/den;
    
    b = (a1+a2)/2;
    
    mse_old = mse_new;
    mse_new = integral(mse_func(a1, a2, b, mean_vec, sigma_vec, alpha_vec), -Inf, Inf);
    iter = iter + 1;
end

partition = b;
codebook = [a1, a2];
p = integral(pdf_func(mean_vec, sigma_vec, alpha_vec), -Inf, b);
output_dist = [1-p p];


end


function z = pdf_channel_output(y, mean_vec, sigma_vec, alpha_vec)
    z = 0 ;
    for k = 1:length(mean_vec)
        z = z + alpha_vec(k) .* 1/sqrt(2*pi*sigma_vec(k)^2) .* exp( -(y-mean_vec(k)).^2/(2*sigma_vec(k)^2) );
    end
end

function z = pdf_func(mean_vec, sigma_vec, alpha_vec)
    z = @(y) pdf_channel_output(y, mean_vec, sigma_vec, alpha_vec);
end

function z = mean_func(mean_vec, sigma_vec, alpha_vec)
    z = @(y) y .* pdf_channel_output(y, mean_vec, sigma_vec, alpha_vec);
end

function z = mse_func(a1, a2, b, mean_vec, sigma_vec, alpha_vec)
    z = @(y) ((y - a1).^2) .* pdf_channel_output(y, mean_vec, sigma_vec, alpha_vec) .* (y < b) + ...
        ((y - a2).^2) .* pdf_channel_output(y, mean_vec, sigma_vec, alpha_vec) .* (y >= b);
end

