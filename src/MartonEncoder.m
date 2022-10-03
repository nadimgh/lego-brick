function [x_out, v_out] = MartonEncoder(u_in, input_dist_enc, H_struct, scl_flag, list_size)
%Encoder for Marton coding scheme
%  u_in: 1xR cell, each containing a vector of information bits to each user
%  input_dist_enc: 1x(2^R) vector representing the input ditribution p(x_1^R) used by encoders
%  H_struct: 1xR cell of structs for the 2R parity check matrices to use (each cell includes both H1 and H2)
%  scl_flag: flag representing whether SCL decoding is used with the polar code
%  list_size: list size to use when SCL decoding is used
%  x_out: 1xR cell of vectors for channel input
%  v_out: 1xR cell, containing vectors of common randomness between the encoder and the decoder

n = size(H_struct{1}.H1, 2);
R = length(u_in);
x_out = cell(1,R);
v_out = cell(1,R);

input_dist = sum(reshape(input_dist_enc, 2, []));
[x_out{1}, v_out{1}] = AsymmetricChannelEncoder(u_in{1}, input_dist, H_struct{1}, scl_flag, list_size);

for i_user = 2:R
    s = zeros(n, i_user-1, 'int8');
    for i_s = 1:i_user-1
        s(:,i_s) = x_out{i_s};
    end
    if i_user < R
        input_dist = sum(reshape(input_dist_enc, 2^(R-i_user), []));
    else
        input_dist = input_dist_enc;
    end
    [x_out{i_user}, v_out{i_user}] = GelfandPinskerEncoder(u_in{i_user}, s, input_dist, H_struct{i_user}, scl_flag, list_size);
end



end

