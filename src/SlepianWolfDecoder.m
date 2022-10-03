function x_out = SlepianWolfDecoder(index, y_in, input_dist, H_struct, subcode_flag, gauss_flag, gauss_params, scl_flag, list_size)
%Slepian-Wolf decoder
%  index: nx1 vector for the shared index from encoder
%  y_in: matrix for the side information (n x num_side_info if gauss_flag == 0 and n x 1 otherwise) 
%  input_dist: ditribution p(y,x) (side information first) if gauss_flag == 0 (not used if gauss_flag == 1)
%  H_struct: struct for parity check matrix used
%  subcode_flag: flag to use either H1 or H2 in H_struct
%  gauss_flag: flag to indicate if p(y|x) is a gaussian mixture or not (if not, then X and Y are taken to be binary)
%  gauss_params: struct that includes (p_x, mean0, sigma0, alpha0, mean1, sigma1, alpha1) for the gaussian mixture parameters given x=0,1 (used only if gauss_flag=1)
%  scl_flag: flag representing whether SCL decoding is used with the polar code
%  list_size: list size to use when SCL decoding is used
%  x_out: reconstruction of source sequence

if ( any(input_dist < 0) || any(input_dist > 1) || abs(sum(input_dist) - 1) > 1e-5 )
    error('Input_dist should be a probability distribution.');
end

if ( length(input_dist) ~= 2^(size(y_in,2)+1) && gauss_flag == 0)
    error('Size of input_dist should be consistent with the number of side information sequences.');
end
input_dist = input_dist(:);

n = length(index);

if subcode_flag == 1
    k = n-size(H_struct.H2, 1);
else
    k = n-size(H_struct.H1, 1); 
end

if subcode_flag == 0
    H = H_struct.H1;
    H_inv = H_struct.H1_inv;
    idx_inv = H_struct.idx1_inv;
else
    H = H_struct.H2;
    H_inv = H_struct.H2_inv;
    idx_inv = H_struct.idx2_inv;
end

v = randi([0 1], n, 1, 'int8');
synV = zeros(n, 1, 'int8');
synV(idx_inv) = mod(H_inv*H*double(v), 2);
r = mod(index + v + synV, 2);

if gauss_flag == 0
	num_side_info = size(y_in, 2);
    if (num_side_info ~= 0)
        y_idx = cast(bi2de(y_in, 'left-msb') + 1, 'int8');
        idx0 = 2*y_idx + r - 1;         % q(y1,y2,v|x) = input_dist(y1,y2,x \xor v)
        idx1 = 2*y_idx - r;
        LLRin = log(input_dist(idx0)./input_dist(idx1));
    else
        alpha = input_dist(2);
        LLRin_alpha = log((1-alpha)/alpha);
        LLRin = zeros(n,1);
        LLRin(r == 0) = LLRin_alpha;
        LLRin(r == 1) = -LLRin_alpha;
    end
else
    mean0 = gauss_params.mean0; sigma0 = gauss_params.sigma0; alpha0 = gauss_params.alpha0;
    mean1 = gauss_params.mean1; sigma1 = gauss_params.sigma1; alpha1 = gauss_params.alpha1;
    p_x = gauss_params.p_x;
    num_mixtures = length(mean0);
    idx0_r = (r == 0);
    idx1_r = (r == 1);
    summ0 = zeros(n,1);
    summ1 = zeros(n,1);
    for i = 1:num_mixtures
        summ0(idx0_r) = summ0(idx0_r) + p_x(1)*alpha0(i)*normpdf(y_in(idx0_r), mean0(i), sigma0(i));
        summ0(idx1_r) = summ0(idx1_r) + p_x(2)*alpha1(i)*normpdf(y_in(idx1_r), mean1(i), sigma1(i));
        summ1(idx0_r) = summ1(idx0_r) + p_x(2)*alpha1(i)*normpdf(y_in(idx0_r), mean1(i), sigma1(i));
        summ1(idx1_r) = summ1(idx1_r) + p_x(1)*alpha0(i)*normpdf(y_in(idx1_r), mean0(i), sigma0(i));
    end
    LLRin = log(summ0./summ1);
end

if subcode_flag
    info_set = H_struct.info_set_2;
else
    info_set = H_struct.info_set_1;
end

if scl_flag == 0 || list_size == 1
    [~, x_est, ~] = polarDecodeSSC(LLRin', info_set);
    x_est = cast(x_est', 'int8');
else
    obj = PolarCode_SCL(n, k, info_set, 0, []);
    [~, x_est, ~] = obj.decode_scl_llr(LLRin, list_size);
    x_est = cast(x_est, 'int8');
end

x_out = mod(x_est + r, 2);

end

