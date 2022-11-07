function [uhat, xhat, yhat] = UplinkCranDecoderCP(index_in, v_user, v_relay, source_decoder1_dist, source_decoder2_dist, channel_decoder1_dist, channel_decoder2_dist, H_struct, scl_flag, list_size)
%Decoder at the central processor for the R-user, T-relay uplink C-RAN problem
%  index_in: 1xT cell of vectors of indices inputted from the relays
%  v_user: 1xR cell, containing vectors of common randomness between the asymmetric channel encoders and decoders
%  v_relay: 1xT cell, containing vectors of common randomness between the source encoders and decoders
%  source_decoder1_dist: marginal ditribution p(y1hat) used by the lossy source decoder
%  source_decoder2_dist: joint ditribution p(y1hat,y2hat) (in this order) used by the Wyner-Ziv decoder
%  channel_decoder1_dist: joint ditribution p(y1hat,y2hat,x1) (in this order) used by the first asymmetric channel decoder
%  channel_decoder2_dist: joint ditribution p(x1,y1hat,y2hat,x2) (in this order) used by the second asymmetric channel decoder
%  H_struct: 1x(R+T) cell of structs for the 2(R+T) parity check matrices to use (each cell includes both H1 and H2)
%  scl_flag: flag representing whether SCL decoding is used with the polar code
%  list_size: list size to use when SCL decoding is used
%  uhat: 1xR cell of vectors for estimated information bits
%  xhat: 1xR cell of vectors for estimated transmitted sequences
%  yhat: 1xT cell of vectors for quantization codewords

T = length(index_in);
R = length(v_user);

uhat = cell(1,R);
xhat = cell(1,R);
yhat = cell(1,T);

yhat{1} = LossySourceDecoder(index_in{1}, v_relay{1}, source_decoder1_dist, H_struct{3}, scl_flag, list_size);
yhat{2} = WynerZivDecoder(index_in{2}, v_relay{2}, yhat{1}, source_decoder2_dist, H_struct{4}, scl_flag, list_size);

y_dec_1 = [yhat{1} yhat{2}];
[uhat{1}, xhat{1}] = AsymmetricChannelDecoder(y_dec_1, v_user{1}, channel_decoder1_dist, H_struct{1}, 0, [], scl_flag, list_size);

y_dec_2 = [xhat{1} yhat{1} yhat{2}];
[uhat{2}, xhat{2}] = AsymmetricChannelDecoder(y_dec_2, v_user{2}, channel_decoder2_dist, H_struct{2}, 0, [], scl_flag, list_size);

end

