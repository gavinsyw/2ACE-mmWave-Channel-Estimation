function variance_with_K_singular_values(H)

[n, ~, ~] = size(H);
cumulative_lam = [];
err = zeros(n, 36);
for i = 1:n
    H_sub = squeeze(H(i,:,:));
    [U, S, V] = svd(H_sub);
    lam = diag(S);
    for j = 1:36
        lam2 = lam;
        lam2(j+1:36) = 0;
        S2 = diag(lam2);
        rec_H = U * S2 * V';
        e = norm(H_sub-rec_H, "fro") / norm(H_sub, "fro");
        err(i, j) = e;
    end
end

plot(1:36, mean(err), "LineWidth", 4);
xlabel("K", "FontSize", 18, "FontName", "Times");
ylabel("Normalized Error", "FontSize", 18, "FontName", "Times");
hold on;
    
end