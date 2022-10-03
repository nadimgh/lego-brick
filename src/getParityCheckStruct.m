function H_cell = getParityCheckStruct(n, rate, bhatta)
%Returns a struct for each entry in (n, rate)
%  n: blocklength
%  rate: (num_codes x 2) vector of rates
%  bhatta: (num_codes x 2) cell of estimated Bhattachatyya parameters of the synthetic channels from the polar transformation (only used when use_polar = 1)
%  H_cell: (1 x num_codes) cell of structs, each containing: (H1, G1, idx1_inv, H1_inv) for first code rate, and (H2, G2, idx2_inv, H2_inv, idx_subcode, idx_subcode_c) 
%  for second code rate, such that G2 is a submatrix of G1 (row indices are stored in idx_subcode), and H1 is a submatrix of H2 (H2 = [H1; Q])

num_codes = size(rate,1);
H_cell = cell(1,num_codes);

m = log2(n);
F = [1 0; 1 1];
G = F;
for i = 1:m-1
    G = kron(G,F);
end
for i = 1:num_codes
    k1 = n*rate(i,1); k2 = n*rate(i,2);
    [~, idx_sorted_1] = sort(bhatta{i,1});

    A_1 = sort(idx_sorted_1(1:k1));
    A_1_c = setdiff(1:n, A_1);
    G1 = G(A_1, :);
    H1 = transpose(G(:, A_1_c));
    idx1_inv = A_1_c;
    H1_inv = getParityCheckSubmatrixInversePolar(H1(:,A_1_c));
    info_set_1 = zeros(1,n);
    info_set_1(A_1) = 1;

    if k2 > 0
        [~, idx_sorted_2] = sort(bhatta{i,2});
        [~, loc] = ismember(idx_sorted_2, A_1);
        idx_sub = loc(find(loc, k2));
        A_2 = sort(A_1(idx_sub));
        A_2_c = setdiff(1:n, A_2);
        A_diff = setdiff(A_1, A_2);
        G2 = G(A_2, :);
        H2 = [H1; transpose(G(:, A_diff))];
        idx2_inv = A_2_c;
        [~, idx_subcode] = ismember(A_2, A_1);
        idx_subcode = sort(idx_subcode);
        idx_subcode_c = setdiff(1:k1, idx_subcode);
        H2_original = transpose(G(:,A_2_c));
        H2_inv_original = getParityCheckSubmatrixInversePolar(H2_original(:,A_2_c));
        [~, permutation] = ismember([A_1_c A_diff], A_2_c);
        H2_inv = H2_inv_original(:, permutation);
        info_set_2 = zeros(1,n);
        info_set_2(A_2) = 1;
    else % k2 = 0
        info_set_2 = zeros(1,n);
        idx2_inv = 1:n;
        G2 = G([], :);
        H2 = [H1; transpose(G(:, A_1))];
        idx_subcode = [];
        idx_subcode_c = 1:k1;
        H2_inv_original = transpose(G);
        [~, permutation] = ismember([A_1_c A_1], 1:n);
        H2_inv = H2_inv_original(:, permutation);
    end

    H_struct = struct('H1', H1, 'G1', G1, 'idx1_inv', idx1_inv, 'H1_inv', H1_inv, 'H2', H2, 'G2', G2, 'idx2_inv', idx2_inv, 'H2_inv', H2_inv, 'idx_subcode', idx_subcode, ...
        'idx_subcode_c', idx_subcode_c, 'info_set_1', info_set_1, 'info_set_2', info_set_2);
    H_cell{i} = H_struct;
end

end

