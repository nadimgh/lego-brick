function H_inv = getParityCheckSubmatrixInversePolar(H)
%get the inverse of a lower-triangular binary square matrix H

r = size(H,1);
H_inv = zeros(r);
for i = r:-1:1
    e = zeros(1,r); e(i) = 1;
    H_inv(i,:) = mod(H(i,i+1:r)*H_inv(i+1:r,:) + e, 2);
end

end

