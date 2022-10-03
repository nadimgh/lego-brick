function [x_out, v1] = AsymmetricChannelEncoder(u_in, input_dist_enc, H_struct, scl_flag, list_size)
%Asymmetric channel encoder
%  u_in: (k1-k2)x1 vector for the information bits
%  input_dist_enc: 1x2 vector for the input distribution used by the encoder
%  H_struct: struct for parity check matrix to use (includes both H1 and H2)
%  scl_flag: flag representing whether SCL decoding is used with the polar code
%  list_size: list size to use when SCL decoding is used
%  x_out: n x 1 vector for channel input
%  v1: (n-k1)x1 vector for common randomness between the encoder and the decoder

n = size(H_struct.H1, 2);
k1 = n-size(H_struct.H1, 1);
k2 = n-size(H_struct.H2, 1);

if (length(u_in) ~= k1-k2)
    error('Length of information bits should be consistent with H_struct.');
end

v1 = randi([0 1], n-k1, 1, 'int8');
H2_inv = H_struct.H2_inv;
idx2_inv = H_struct.idx2_inv;

index_slepian = zeros(n, 1, 'int8');
index_slepian(idx2_inv) = mod(H2_inv*double([v1; u_in]), 2);

x_out = SlepianWolfDecoder(index_slepian, [], input_dist_enc, H_struct, 1, 0, [], scl_flag, list_size);

end

