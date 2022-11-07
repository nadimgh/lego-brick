function [v_out, x_out] = UplinkCranEncoderUser(u_in, encoder1_dist, encoder2_dist, H_struct, scl_flag, list_size)
%Encoder at the users for the R-user, T-relay uplink C-RAN problem
%  u_in: 1xR cell, each containing a vector of information bits to each user
%  encoder1_dist: joint ditribution p(y1hat,y2hat,x1) (in this order) used by the first asymmetric channel encoder
%  encoder2_dist: joint ditribution p(x1,y1hat,y2hat,x2) (in this order) used by the second asymmetric channel encoder
%  H_struct: 1x(R+T) cell of structs for the 2(R+T) parity check matrices to use (each cell includes both H1 and H2)
%  scl_flag: flag representing whether SCL decoding is used with the polar code
%  list_size: list size to use when SCL decoding is used
%  v_out: 1xR cell, containing vectors of common randomness between the asymmetric channel encoders and decoders
%  x_out: 1xR cell of vectors of channel inputs

R = length(u_in);

x_out = cell(1,R);
v_out = cell(1,R);

H_struct_user = cell(1,R);
for i = 1:R
    H_struct_user{i} = H_struct{i};
end
[x_out{1}, v_out{1}] = AsymmetricChannelEncoder(u_in{1}, encoder1_dist, H_struct{1}, scl_flag, list_size);
[x_out{2}, v_out{2}] = AsymmetricChannelEncoder(u_in{2}, encoder2_dist, H_struct{2}, scl_flag, list_size);

end

