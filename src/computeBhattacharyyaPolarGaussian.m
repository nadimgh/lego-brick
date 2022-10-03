function bhatta = computeBhattacharyyaPolarGaussian(n, input_dist, gauss_params, symmetrized_flag)
%Computes estimated Bhattacharyya parameters of the synthetic channels from
%the polar transformation for a Gaussian channel
%  input_dist: channel input distribution in the form [1-alpha alpha] (used only if symmetrized flag=1)
%  gauss_params: struct that includes (mean0, sigma0, alpha0, mean1, sigma1, alpha1) for the Gaussian channel parameters given x=0,1
%  symmetrized_flag: flag indicating if Bhattacharyya parameters corresponding to the symmetrized channel of the input-output joint distribution are needed

if (n == 1)
    mean0 = gauss_params.mean0; sigma0 = gauss_params.sigma0; alpha0 = gauss_params.alpha0;
    mean1 = gauss_params.mean1; sigma1 = gauss_params.sigma1; alpha1 = gauss_params.alpha1;
    low_val = min([mean0 mean1])-5*max([sigma0 sigma1]);
    upper_val = max([mean0 mean1])+5*max([sigma0 sigma1]);
    bhatta = integral(bhattacharyya_function(mean0, sigma0, alpha0, mean1, sigma1, alpha1), low_val, upper_val);
    if symmetrized_flag == 1
        bhatta = 2*sqrt(input_dist(1)*input_dist(2))*bhatta;
    end
else
    z = computeBhattacharyyaPolarGaussian(n/2, input_dist, gauss_params, symmetrized_flag);
    bhatta = [2*z-z.^2; z.^2];
    bhatta = bhatta(:)';
end

end


function z = pdf_channel_output(y, mean, sigma, alpha)
    z = 0 ;
    for k = 1:length(mean)
        z = z + alpha(k) .* 1/sqrt(2*pi*sigma(k)^2) .* exp( -(y-mean(k)).^2/(2*sigma(k)^2) );
    end
end

function z = bhattacharyya_function(mean0, sigma0, alpha0, mean1, sigma1, alpha1)
    z = @(y) sqrt( pdf_channel_output(y, mean0, sigma0, alpha0) .* pdf_channel_output(y, mean1, sigma1, alpha1) );
end

