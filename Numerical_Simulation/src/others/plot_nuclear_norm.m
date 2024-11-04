load("Synthetic_result_Vs_M_L1_SNR_Inf_TX12_RX12.mat");
% H = load(["../Dataset/H_livingroom_180Rx.mat"]);
% H = H.H;

H_Estimated = Simulation_result_Vs_M.H_Estimated;
Mean_Evaluation = Simulation_result_Vs_M.Mean_Evaluation;
H = Simulation_result_Vs_M.H;

H_PLGAMP = squeeze(H_Estimated(1,:,5,:,:));
H_PLOMP = squeeze(H_Estimated(1,:,6,:,:));
H_ADMM = squeeze(H_Estimated(1,:,12,:,:));
H_PhaseLift = squeeze(H_Estimated(1,:,1,:,:));

H_gt = H.H_Matrix;
% H_gt = squeeze(H(1,:,:));
% H_gt = H_gt / norm(squeeze(H_gt), "fro");
% M = 6:2:23;

norm_PLGAMP = zeros(numel(M), 1);
norm_PLOMP = zeros(numel(M), 1);
norm_ADMM = zeros(numel(M), 1);
norm_PhaseLift = zeros(numel(M), 1);

for M_idx = 1:numel(M)
    H_tmp = squeeze(H_PLOMP(M_idx, :, :));
    [U,S,V] = svd(H_tmp);
    norm_PLOMP(M_idx) = sum(abs(diag(S)));
    H_tmp = squeeze(H_PLGAMP(M_idx, :, :));
    [U,S,V] = svd(H_tmp);
    norm_PLGAMP(M_idx) = sum(abs(diag(S)));
    H_tmp = squeeze(H_ADMM(M_idx, :, :));
    [U,S,V] = svd(H_tmp);
    norm_ADMM(M_idx) = sum(abs(diag(S)));
    H_tmp = squeeze(H_PhaseLift(M_idx,:,:));
    [U,S,V] = svd(H_tmp);
    norm_PhaseLift(M_idx) = sum(abs(diag(S)));
end

[U,S,V] = svd(H_gt);
norm_gt = sum(abs(diag(S)));

figure();
plot(M.^2, norm_PLOMP, "LineWidth", 4);
hold on;
plot(M.^2, norm_PLGAMP, "LineWidth", 4);
% hold on;
% plot(M.^2, norm_ADMM, "LineWidth", 4);
hold on;
plot(M.^2, norm_PhaseLift, "LineWidth", 4);
hold on;
plot(M.^2, ones(numel(M), 1)*norm_gt, '--', "LineWidth", 4);

grid on;
xlabel("T", "FontSize", 20, "Interpreter", 'latex');
ylabel("Nuclear Norm", "FontSize", 20, "Interpreter", 'latex');
legend("PLOMP", "PLGAMP", "PhaseLift", "Ground Truth", "FontSize", 20);
ax = gca;
ax.FontSize = 20;
