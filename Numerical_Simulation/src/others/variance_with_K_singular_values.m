function err = variance_with_K_singular_values(H, plot_flag)

[n, TX, ~] = size(H);
err = zeros(n, TX);
for i = 1:n
    H_sub = squeeze(H(i,:,:));
    [U, S, V] = svd(H_sub);
    lam = diag(S);
    for j = 1:TX-1
        lam2 = lam;
        lam2(j+1:TX) = 0;
%         S2 = diag(lam2);
        e = sum(abs(lam2).^2) / sum(abs(lam).^2);
        err(i, j) = e;
    end
end

% standard error assumption
f0 = 0.8;
f1 = 0.9;
f2 = 0.95;
f3 = 0.995;
std_err_A1 = zeros(1, TX);
std_err_A1(1.5*ceil(sqrt(TX)):TX) = 0.95;
std_err_A2 = zeros(1, TX);
std_err_A2(ceil(0.7*sqrt(TX)):TX) = f0;
std_err_A2(ceil(1*sqrt(TX)):TX) = f1;
std_err_A2(ceil(1.5*sqrt(TX)):TX) = f2;
std_err_A2(ceil(3*sqrt(TX)):TX) = f3;

if plot_flag
    figure();
%     plot(log(1:TX), pow2db(1-transpose(err)), "LineWidth", 4);
    loglog(1:TX, transpose(err), "LineWidth", 4);
    hold on;
%     plot(log(1:TX), pow2db(std_err_A1), "g--", "LineWidth", 4);
    loglog(1:TX, std_err_A1, "g--", "LineWidth", 4);
    hold on;
%     plot(log(1:TX), pow2db(std_err_A2), "y--", "LineWidth", 4);
    loglog(1:TX, std_err_A2, "r--", "LineWidth", 4);
    grid on;
    xticks([1 4 9 16 25 36]);
    xlabel("K", "FontSize", 20, "Interpreter", 'latex');
    ylabel("Captured Energy", "FontSize", 20, "Interpreter", 'latex');
%     legend("Trace 1", "Trace 2", "Trace 3", "Trace 4", "Trace 5", "A1 Profile", "A2 Profile", "FontSize", 20);
    ax = gca;
    ax.FontSize = 20;
end
    
end