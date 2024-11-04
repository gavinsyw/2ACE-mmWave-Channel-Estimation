function n = nuclear_norm(H)

[n, ~, ~] = size(H);
nc_norm = zeros(n);
for i = 1:n
    H_sub = squeeze(H(i,:,:));
    V = svd(H_sub);
    nc_norm = sum(V) ./ norm(H_sub, "fro");
end

n = mean(nc_norm);

end