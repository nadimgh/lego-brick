function [ A ] = get_tuples( n )
%get all binary n-tuples

D = 0:2^n - 1;
B = dec2bin(D);
A = zeros(2^n,n);
for i=1:n
    A(:,i) = str2num(B(:,i));
end
end
