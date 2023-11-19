function CDF_H(H)

[n, ~, ~] = size(H);
cumulative_lam = [];
for i = 1:n
    H_sub = squeeze(H(i,:,:));
    [~,lam] = eig(H_sub);
    lam = abs(diag(lam));
    lam = lam ./ max(lam);
    cumulative_lam = [cumulative_lam; lam];
end

cdfplot(cumulative_lam);
xlabel("Normalized eigenvalue", "FontSize", 18, "FontName", "Times");
ylabel("CDF", "FontSize", 18, "FontName", "Times");

end