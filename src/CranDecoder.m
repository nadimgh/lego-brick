function u_out = CranDecoder(y_in, v_in, H_struct, gauss_params, scl_flag, list_size)
%Decoder for the R-user, T-relay C-RAN problem
%  y_in: 1xR cell of vectors for the channel output
%  v_in: 1x(R+T) cell, containing vectors of common randomness between the encoder of CP, encoders at relays and decoders
%  H_struct: 1x(R+T) cell of structs for the 2(R+T) parity check matrices to use (each cell includes both H1 and H2)
%  gauss_params: 1xR cell of structs that include (mean0, sigma0, alpha0, mean1, sigma1, alpha1) for the gaussian parameters given x_i=0,1 for each i
%  scl_flag: flag representing whether SCL decoding is used with the polar code
%  list_size: list size to use when SCL decoding is used
%  u_out: 1xR cell of vectors for estimated information bits

R = length(y_in);

H_struct_marton = cell(1,R);
for i = 1:R
    H_struct_marton{i} = H_struct{i};
end

v_marton = cell(1,R);
for i = 1:R
    v_marton{i} = v_in{i};
end

u_out = MartonDecoder(y_in, v_marton, H_struct_marton, gauss_params, scl_flag, list_size);

end

