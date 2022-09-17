function bhatta = computeBhattacharyyaPolar(n, target_dist)
%Computes estimated Bhattacharyya parameters of the synthetic channels from
%a symmetrized channel with joint distribution target_dist

if (n == 1)
    target_dist_reshape = reshape(target_dist, 2, []);
    bhatta = 2*sum(sqrt(target_dist_reshape(1,:).*target_dist_reshape(2,:)));
else
    z = computeBhattacharyyaPolar(n/2, target_dist);
    bhatta = [2*z-z.^2; z.^2];
    bhatta = bhatta(:)';
end

return

