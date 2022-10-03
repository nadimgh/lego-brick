function index = SlepianWolfEncoder(x_in, H_struct, subcode_flag)
%Slepian-Wolf encoder

n = length(x_in);
if subcode_flag == 0
    H = H_struct.H1;
    H_inv = H_struct.H1_inv;
    idx_inv = H_struct.idx1_inv;
else
    H = H_struct.H2;
    H_inv = H_struct.H2_inv;
    idx_inv = H_struct.idx2_inv;
end

index = zeros(n, 1, 'int8');
index(idx_inv) = mod(H_inv*H*double(x_in), 2);

end

