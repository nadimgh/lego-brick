function u_out = MartonDecoder(y_in, v_in, H_struct, gauss_params, scl_flag, list_size)
%Decoder for Marton coding scheme
%  y_in: 1xR cell of vectors for the channel output
%  v_in: 1xR cell of vectors for common randomness between the encoder and the decoder
%  H_struct: 1xR cell of structs for the 2R parity check matrices to use (each cell includes both H1 and H2)
%  gauss_params: 1xR cell of structs that include (mean0, sigma0, alpha0, mean1, sigma1, alpha1) for the gaussian parameters given x_i=0,1 for each i
%  scl_flag: flag representing whether SCL decoding is used with the polar code
%  list_size: list size to use when SCL decoding is used
%  u_out: 1xR cell of vectors for estimated information bits

R = length(y_in);
u_out = cell(1,R);

u_out{1} = AsymmetricChannelDecoder(y_in{1}, v_in{1}, 1, H_struct{1}, 1, gauss_params{1}, scl_flag, list_size);

for i_user = 2:R
    u_out{i_user} = GelfandPinskerDecoder(y_in{i_user}, v_in{i_user}, 1, H_struct{i_user}, 1, gauss_params{i_user}, scl_flag, list_size);
end

end

