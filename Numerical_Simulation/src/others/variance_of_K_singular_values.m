function err = variance_of_K_singular_values(H, plot_flag)

[n, TX, ~] = size(H);
err = zeros(n, TX);
for i = 1:n
    H_sub = squeeze(H(i,:,:));
    [U, S, V] = svd(H_sub);
    lam = diag(S);
    for j = 1:TX
        spot_eig = lam(j);
        lam2 = zeros(1, TX);
        lam2(j) = spot_eig;
        S2 = diag(lam2);
        e = abs(spot_eig)^2 / sum(abs(lam).^2);
        err(i, j) = e;
    end
end

if plot_flag
    figure();
%     plot(log(1:TX), pow2db(1-transpose(err)), "LineWidth", 4);
    loglog(1:TX, transpose(err), "LineWidth", 4);
    grid on;
    xticks([1 4 9 16 25 36]);
    xlabel("K", "FontSize", 20, "Interpreter", 'latex');
    ylabel("Captured Energy", "FontSize", 20, "Interpreter", 'latex');
    legend("Trace 1", "Trace 2", "Trace 3", "Trace 4", "Trace 5", "FontSize", 20);
    ax = gca;
    ax.FontSize = 20;
end
    
end