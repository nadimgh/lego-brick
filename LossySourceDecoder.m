function x_out = LossySourceDecoder(index_in, v_in, input_dist_dec, H_struct, scl_flag, list_size)
%Lossy source decoder
%  index_in: (k1-k2)x1 vector for the index inputted to the decoder
%  v_in: (n-k1)x1 vector for common randomness between the encoder and the decoder
%  input_dist_dec: the marginal ditribution p(\hat{x}) used by the decoder (should have the form [1-alpha alpha])
%  H_struct: struct for parity check matrix to use (includes both H1 and H2)
%  scl_flag: flag representing whether SCL decoding is used with the polar code
%  list_size: list size to use when SCL decoding is used
%  x_out: nx1 vector for the reconstruction of the source sequence

n = size(H_struct.H1, 2);

idx2_inv = H_struct.idx2_inv;
H2_inv = H_struct.H2_inv;

index_slepian = zeros(n, 1, 'int8');
index_slepian(idx2_inv) = mod(H2_inv*double([v_in; index_in]), 2);

x_out = SlepianWolfDecoder(index_slepian, [], input_dist_dec, H_struct, 1, 0, [], scl_flag, list_size);

end

