function eig_decay(H)

[n, Nt, ~] = size(H);
lam_dis = zeros(n, Nt);
for i = 1:n
    H_sub = squeeze(H(i,:,:));
    lam = abs(eig(H_sub));
    lam = lam ./ max(lam);
    lam_dis(i,:) = log(lam);
end

boxplot(lam_dis);
xlabel("Eigenvalue index", "FontSize", 18, "FontName", "Times");
ylabel("Relative log value", "FontSize", 18, "FontName", "Times");

end