function checkInputParameters_UplinkCran(n, rate, dist_struct, H, C_vec)
%check if the parameters of the uplink C-RAN coding scheme have the correct form

asym_encoder1_dist = dist_struct.asym_encoder1;
asym_encoder2_dist = dist_struct.asym_encoder2;
asym_decoder1_dist = dist_struct.asym_decoder1;
asym_decoder2_dist = dist_struct.asym_decoder2;
lossy_encoder_dist = dist_struct.lossy_encoder;
lossy_decoder_dist = dist_struct.lossy_decoder;
wz_encoder_dist = dist_struct.wz_encoder;
wz_decoder_dist = dist_struct.wz_decoder;
num_power = size(asym_encoder1_dist,1);
[T, R] = size(H); % R-user, T-relay

for i_p = 1:num_power
    if ( any(asym_encoder1_dist(i_p,:) < 0) || any(asym_encoder1_dist(i_p,:) > 1) || abs(sum(asym_encoder1_dist(i_p,:)) - 1) > 1e-5 )
        error('Each row in dist_struct should be a probability distribution.');
    end
end

for i_p = 1:num_power
    if ( any(asym_encoder2_dist(i_p,:) < 0) || any(asym_encoder2_dist(i_p,:) > 1) || abs(sum(asym_encoder2_dist(i_p,:)) - 1) > 1e-5 )
        error('Each row in dist_struct should be a probability distribution.');
    end
end

for i_p = 1:num_power
    if ( any(asym_decoder1_dist(i_p,:) < 0) || any(asym_decoder1_dist(i_p,:) > 1) || abs(sum(asym_decoder1_dist(i_p,:)) - 1) > 1e-5 )
        error('Each row in dist_struct should be a probability distribution.');
    end
end

for i_p = 1:num_power
    if ( any(asym_decoder2_dist(i_p,:) < 0) || any(asym_decoder2_dist(i_p,:) > 1) || abs(sum(asym_decoder2_dist(i_p,:)) - 1) > 1e-5 )
        error('Each row in dist_struct should be a probability distribution.');
    end
end

for i_p = 1:num_power
    if ( any(lossy_encoder_dist(i_p,:) < 0) || any(lossy_encoder_dist(i_p,:) > 1) || abs(sum(lossy_encoder_dist(i_p,:)) - 1) > 1e-5 )
        error('Each row in dist_struct should be a probability distribution.');
    end
end

for i_p = 1:num_power
    if ( any(lossy_decoder_dist(i_p,:) < 0) || any(lossy_decoder_dist(i_p,:) > 1) || abs(sum(lossy_decoder_dist(i_p,:)) - 1) > 1e-5 )
        error('Each row in dist_struct should be a probability distribution.');
    end
end

for i_p = 1:num_power
    if ( any(wz_encoder_dist(i_p,:) < 0) || any(wz_encoder_dist(i_p,:) > 1) || abs(sum(wz_encoder_dist(i_p,:)) - 1) > 1e-5 )
        error('Each row in dist_struct should be a probability distribution.');
    end
end

for i_p = 1:num_power
    if ( any(wz_decoder_dist(i_p,:) < 0) || any(wz_decoder_dist(i_p,:) > 1) || abs(sum(wz_decoder_dist(i_p,:)) - 1) > 1e-5 )
        error('Each row in dist_struct should be a probability distribution.');
    end
end

if (size(asym_decoder2_dist,2) ~= 2^(R+T))
    error('The number of columns in the target distribution of the second asymmetric channel code should be equal to 2^K, where K is the total number of users and relays in the C-RAN problem.');
end

if ( any(C_vec < 0) || any(C_vec > 1) || length(C_vec) ~= T )
    error('Each entry in C_vec should be between zero and one, and the length of C_vec should be equal to the number of relays.')
end

if (size(rate,1) ~= R+T)
    error('The number of rows in rate should be equal to the number of users and relays in the C-RAN problem.');
end

if any(n*rate - round(n*rate) > 1e-5, 'all')
    error('n*rate should be an integer.');
end

for i_relay = R+1:R+T
    if rate(i_relay,1) - rate(i_relay,2) > C_vec(i_relay-R)
        error('Compression rate should be smaller than the backhaul constraints.');
    end
end



end

