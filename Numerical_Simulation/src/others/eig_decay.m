function eig_decay(H)

[n, Nt, ~] = size(H);
lam_dis = zeros(n, Nt);
for i = 1:n
    H_sub = squeeze(H(i,:,:));
    lam = abs(eig(H_sub));
    lam = lam ./ max(lam);
    lam_dis(i,:) = lam;
end

plot(log(mean(lam_dis)), "LineWidth", 4);
xlabel("Eigenvalue index", "FontSize", 18, "Interpreter", 'latex');
ylabel("Relative log value", "FontSize", 18, "Interpreter", 'latex');

end