function checkInputParameters_Cran(n, rate, target_dist, H, C_vec)
%check if the input parameters to sim_cran_marton have the correct form

num_power = size(target_dist,1);
[R, T] = size(H);

for i_p = 1:num_power
    if ( any(target_dist(i_p,:) < 0) || any(target_dist(i_p,:) > 1) || abs(sum(target_dist(i_p,:)) - 1) > 1e-5 )
        error('Each row in target_dist should be a probability distribution.');
    end
end

if (size(target_dist,2) ~= 2^(R+T))
    error('The number of columns in target_dist should be equal to 2^K, where K is the total number of users and relays in the C-RAN problem.');
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
        error('Compression rate should be smaller than the fronthaul constraints.');
    end
end



end
