function [x_out, v1] = GelfandPinskerEncoder(u_in, s_in, input_dist_enc, H_struct, scl_flag, list_size)
%Asymmetric channel encoder
%  u_in: (k1-k2)x1 vector for the information bits
%  s_in: (n x num_states) vector for the state sequence
%  input_dist_enc: target ditribution p(s,x) (channel state first, i.e., has form [P(S=0,X=0) P(S=0,X=1) P(S=1,X=0) P(S=1,X=1)])
%  H_struct: struct for parity check matrix to use (includes both H1 and H2)
%  scl_flag: flag representing whether SCL decoding is used with the polar code
%  list_size: list size to use when SCL decoding is used
%  x_out: n x 1 vector for channel input
%  v1: (n-k1)x1 vector for common randomness between the encoder and the decoder

n = size(H_struct.H1, 2);
k1 = n-size(H_struct.H1, 1);
k2 = n-size(H_struct.H2, 1);
num_states = size(s_in, 2);

if (length(u_in) ~= k1-k2)
    error('Length of information bits should be consistent with H_struct.');
end

v1 = randi([0 1], n-k1, 1, 'int8');
H2_inv = H_struct.H2_inv;
idx2_inv = H_struct.idx2_inv;

p_s = sum(reshape(input_dist_enc, 2, []));
beta = input_dist_enc(2)/p_s(1);  % P(X=1 | S=0)
gamma = input_dist_enc(3)/p_s(2); % P(X=0 | S=1)

index_slepian = zeros(n, 1, 'int8');
index_slepian(idx2_inv) = mod(H2_inv*double([v1; u_in]), 2);

if beta == gamma && num_states == 1
    x_out = SlepianWolfDecoder(mod(index_slepian + s_in, 2), [], [1-beta beta], H_struct, 1, 0, [], scl_flag, list_size);
    x_out = mod(x_out + s_in, 2);
else
    x_out = SlepianWolfDecoder(index_slepian, s_in, input_dist_enc, H_struct, 1, 0, [], scl_flag, list_size);
end

end

