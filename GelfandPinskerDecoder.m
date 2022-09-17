function u_out = GelfandPinskerDecoder(y_in, v_in, input_dist_dec, H_struct, gauss_flag, gauss_params, scl_flag, list_size)
%Gelfand-Pinsker decoder
%  y_in: n x 1 vector for the channel output
%  v_in: (n-k1)x1 vector for common randomness between the encoder and the decoder
%  input_dist_dec: if gauss_flag == 0, this gives ditribution p(s,y,x) (state distribution first, then channel output, then channel input). Not used if gauss_flag == 1
%  H_struct: struct for parity check matrix to use (includes both H1 and H2)
%  gauss_flag: flag to indicate if p(y|x) is a gaussian mixture or not (if not, then X and Y are taken to be binary)
%  gauss_params: struct that includes (p_x, mean0, sigma0, alpha0, mean1, sigma1, alpha1) for the gaussian mixture parameters given x=0,1 (used only if gauss_flag=1)
%  scl_flag: flag representing whether SCL decoding is used with the polar code
%  list_size: list size to use when SCL decoding is used
%  u_out: k1 x 1 vector for the estimated information bits

n = size(y_in, 1);
k1 = n-size(H_struct.H1, 1);

H2 = H_struct.H2;
H1_inv = H_struct.H1_inv;
idx1_inv = H_struct.idx1_inv;

v = zeros(n, 1, 'int8');
v(idx1_inv) = mod(H1_inv*double(v_in), 2);

if gauss_flag == 0
    input_dist_slepian = sum(reshape(input_dist_dec, 4, []), 2);
else
    input_dist_slepian = 1;
end

xhat = SlepianWolfDecoder(v, y_in, input_dist_slepian, H_struct, 0, gauss_flag, gauss_params, scl_flag, list_size);
u_est = mod(H2*double(xhat), 2);
u_out = u_est(n-k1+1:end);

end