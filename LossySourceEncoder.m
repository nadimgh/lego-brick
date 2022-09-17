function [index_out, v_out, x_out] = LossySourceEncoder(x_in, input_dist_enc, H_struct, scl_flag, list_size)
%Lossy source ncoder
%  x_in: nx1 vector for the source sequence
%  input_dist_enc: 1x4 vector for the joint distribution p(x,\hat{x}) used by the encoder
%  H_struct: struct for parity check matrix to use (includes both H1 and H2)
%  scl_flag: flag representing whether SCL decoding is used with the polar code
%  list_size: list size to use when SCL decoding is used
%  index_out: nx1 vector for the index shared to the decoder
%  v_out: (n-k1)x1 vector for common randomness between the encoder and the decoder
%  x_out: nx1 vector for the sequence compressed by the encoder and the decoder

n = size(x_in, 1);

k1 = n-size(H_struct.H1, 1);
idx1_inv = H_struct.idx1_inv;
H1_inv = H_struct.H1_inv;
Q = H_struct.H2(n-k1+1:end,:);

v_out = randi([0 1], n-k1, 1, 'int8');
v_slepian = zeros(n, 1, 'int8');
v_slepian(idx1_inv) = mod(H1_inv*double(v_out), 2);

x_out = SlepianWolfDecoder(v_slepian, x_in, input_dist_enc, H_struct, 0, 0, [], scl_flag, list_size);
index_out = mod(Q*double(x_out), 2);

end

