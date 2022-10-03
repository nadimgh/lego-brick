function [index_out, v_out, x_out] = CranEncoderCP(u_in, target_dist, H_struct, scl_flag, list_size)
%Encoder at the CP for the R-user, T-relay C-RAN problem
%  u_in: 1xR cell, each containing a vector of information bits to each user
%  target_dist: 1x2^(R+T) vector for the target ditribution p(u_1^R,x_1^T) used by the Marton encoder and the compressor
%  H_struct: 1x(R+T) cell of structs for the 2(R+T) parity check matrices to use (each cell includes both H1 and H2)
%  scl_flag: flag representing whether SCL decoding is used with the polar code
%  list_size: list size to use when SCL decoding is used
%  index_out: 1xT cell of vectors of indices to the relays
%  v_out: 1x(R+T) cell, containing vectors of common randomness between the encoder of CP and encoders at relays
%  x_out: 1x(R+T) cell of vectors of channel inputs generated at the CP

n = size(H_struct{1}.H1, 2);
R = length(u_in);
T = length(H_struct) - R;

target_dist_marton = sum(reshape(target_dist, 2^R, []));

H_struct_marton = cell(1,R);
for i = 1:R
    H_struct_marton{i} = H_struct{i};
end

[x_marton, v_marton] = MartonEncoder(u_in, target_dist_marton, H_struct_marton, scl_flag, list_size);

x_out = cell(1,R+T);
for i = 1:R
    x_out{i} = x_marton{i};
end

v_out = cell(1,R+T);
for i = 1:R
    v_out{i} = v_marton{i};
end

index_out = cell(1,T);
for i_user = R+1:R+T
    k1 = n-size(H_struct{i_user}.H1, 1);
    idx1_inv = H_struct{i_user}.idx1_inv;
    H1_inv = H_struct{i_user}.H1_inv;
    Q = H_struct{i_user}.H2(n-k1+1:end,:);
    
    v_out{i_user} = randi([0 1], n-k1, 1, 'int8');
    v_slepian = zeros(n, 1, 'int8');
    v_slepian(idx1_inv) = mod(H1_inv*double(v_out{i_user}), 2);
    
    s_in = zeros(n, i_user-1);
    for i_s = 1:i_user-1
        s_in(:,i_s) = x_out{i_s};
    end
    
    if i_user < R+T
        target_dist_slepian = sum(reshape(target_dist, 2^(R+T-i_user), []));
    else
        target_dist_slepian = target_dist;
    end
    
    x_out{i_user} = SlepianWolfDecoder(v_slepian, s_in, target_dist_slepian, H_struct{i_user}, 0, 0, [], scl_flag, list_size);
    index_out{i_user-R} = mod(Q*double(x_out{i_user}), 2);
end

end

