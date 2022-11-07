function [v_out, index_out] = UplinkCranEncoderRelay(y_in, encoder1_dist, encoder2_dist, H_struct, scl_flag, list_size)
%Encoder at the relays for the R-user, T-relay uplink C-RAN problem
%  y_in: 1xT cell, each containing a vector of quantized bits at each relay
%  encoder1_dist: joint ditribution p(y1,y1hat) (in this order) used by the lossy source encoder
%  encoder2_dist: joint ditribution p(y2,y2hat) (in this order) used by the Wyner-Ziv encoder
%  H_struct: 1x(R+T) cell of structs for the 2(R+T) parity check matrices to use (each cell includes both H1 and H2)
%  scl_flag: flag representing whether SCL decoding is used with the polar code
%  list_size: list size to use when SCL decoding is used
%  v_out: 1xT cell, containing vectors of common randomness between the encoders and decoders
%  index_out: 1xT cell of vectors of indices to the central processor

T = length(y_in);
R = length(H_struct) - T;

v_out = cell(1,T);
index_out = cell(1,T);

H_struct_relay = cell(1,T);
for i = R+1:R+T
    H_struct_relay{i-R} = H_struct{i};
end
[index_out{1}, v_out{1}, ~] = LossySourceEncoder(y_in{1}, encoder1_dist, H_struct_relay{1}, scl_flag, list_size);
[index_out{2}, v_out{2}, ~] = WynerZivEncoder(y_in{2}, encoder2_dist, H_struct_relay{2}, scl_flag, list_size);

end

