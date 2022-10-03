function x_out = CranEncoderRelay(index_in, v_in, target_dist_relay, H_struct, scl_flag, list_size)
%Encoder at the relays for the R-user, T-relay C-RAN problem
%  index_in: 1xT cell of vectors of indices inputted to the relays
%  v_in: 1x(R+T) cell, containing vectors of common randomness between the encoder of CP, encoders at relays and decoders
%  target_dist_relay: 1xT cell, where i-th cell is the marginal ditribution p(x_{i+R}) (should have the form [1-alpha alpha])
%  H_struct: 1x(R+T) cell of structs for the 2(R+T) parity check matrices to use (each cell includes both H1 and H2)
%  scl_flag: flag representing whether SCL decoding is used with the polar code
%  list_size: list size to use when SCL decoding is used
%  x_out: 1xT cell of vectors of channel inputs

n = size(H_struct{1}.H1, 2);
T = length(index_in);
R = length(H_struct) - T;

x_out = cell(1,T);
for i_user = R+1:R+T
    idx2_inv = H_struct{i_user}.idx2_inv;
    H2_inv = H_struct{i_user}.H2_inv;
    
    index_slepian = zeros(n, 1, 'int8');
    index_slepian(idx2_inv) = mod(H2_inv*double([v_in{i_user}; index_in{i_user-R}]), 2);
    
    x_out{i_user-R} = SlepianWolfDecoder(index_slepian, [], target_dist_relay{i_user-R}, H_struct{i_user}, 1, 0, [], scl_flag, list_size);
end

end

