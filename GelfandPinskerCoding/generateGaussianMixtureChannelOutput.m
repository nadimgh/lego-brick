function y = generateGaussianMixtureChannelOutput(x, s, gauss_params)
%generate output of a Gaussian mixture channel with inputs x and s, 
%where: gauss_params: struct that includes (mean00, sigma00, alpha00, mean01, sigma01, alpha01, mean10, sigma10, alpha10, mean11, sigma11, alpha11) for the gaussian mixture parameters given s,x=0,1 (in this order)

n = length(x);
idx00 = (s==0 & x==0); idx01 = (s==0 & x==1); idx10 = (s==1 & x==0); idx11 = (s==1 & x==1);
num00 = sum(idx00); num01 = sum(idx01); num10 = sum(idx10); num11 = sum(idx11);

idx_smpl = zeros(n, 1);
idx_smpl(idx00) = randsmpl(gauss_params.alpha00, num00, 1, 'int16');
idx_smpl(idx01) = randsmpl(gauss_params.alpha01, num01, 1, 'int16');
idx_smpl(idx10) = randsmpl(gauss_params.alpha10, num10, 1, 'int16');
idx_smpl(idx11) = randsmpl(gauss_params.alpha11, num11, 1, 'int16');

y = zeros(size(x)); 
y(idx00) = gauss_params.mean00(idx_smpl(idx00)) + gauss_params.sigma00(idx_smpl(idx00)).*randn(num00,1);
y(idx01) = gauss_params.mean01(idx_smpl(idx01)) + gauss_params.sigma01(idx_smpl(idx01)).*randn(num01,1);
y(idx10) = gauss_params.mean10(idx_smpl(idx10)) + gauss_params.sigma10(idx_smpl(idx10)).*randn(num10,1);
y(idx11) = gauss_params.mean11(idx_smpl(idx11)) + gauss_params.sigma11(idx_smpl(idx11)).*randn(num11,1);

end

