H_Estimated = Simulation_result_Vs_M.H_Estimated;
Mean_Evaluation = Simulation_result_Vs_M.Mean_Evaluation;
H = Simulation_result_Vs_M.H;

mse_PLGAMP = Mean_Evaluation(1,:,4,6);
mse_ADMM = Mean_Evaluation(1,:,4,12);
mse_PhaseLift = Mean_Evaluation(1,:,4,1);

H_PLGAMP = squeeze(H_Estimated(2,:,6,:,:));
H_ADMM = squeeze(H_Estimated(2,:,12,:,:));
H_PhaseLift = squeeze(H_Estimated(2,:,1,:,:));

H_gt = H.H_Matrix;
M = 6:2:24;

norm_PLGAMP = zeros(numel(M), 1);
norm_ADMM = zeros(numel(M), 1);
norm_PhaseLift = zeros(numel(M), 1);

H_squeezed = zeros(4, 12, 12);
H_squeezed(1,:,:) = squeeze(H_gt);
H_squeezed(2,:,:) = squeeze(H_PLGAMP(5,:,:));
H_squeezed(3,:,:) = squeeze(H_ADMM(5,:,:));
H_squeezed(4,:,:) = squeeze(H_PhaseLift(5,:,:));
plot_flag = 0;

var_gt = variance_with_K_singular_values(H_squeezed, plot_flag);
var_gt = var_gt(1,:);
var_PLGAMP = variance_with_K_singular_values(H_PLGAMP, plot_flag);
var_ADMM = variance_with_K_singular_values(H_ADMM, plot_flag);
var_PhaseLift = variance_with_K_singular_values(H_PhaseLift, plot_flag);

err_PLGAMP = zeros(1, 10);
err_ADMM = zeros(1, 10);
err_PhaseLift = zeros(1, 10);
for i = 1:10
    err_PLGAMP(i) = norm(var_PLGAMP(i,:)-var_gt) / norm(var_gt);
    err_ADMM(i) = norm(var_ADMM(i,:)-var_gt) / norm(var_gt);
    err_PhaseLift(i) = norm(var_PhaseLift(i,:)-var_gt) / norm(var_gt);
end

figure();
scatter(err_PLGAMP, mse_PLGAMP);
hold on;
scatter(err_ADMM, mse_ADMM);
hold on;
scatter(err_PhaseLift, mse_PhaseLift);
grid on;
xlabel("Deviation from Power Law", "FontSize", 20, "Interpreter", 'latex');
ylabel("MSE", "FontSize", 20, "Interpreter", 'latex');
legend("PLGAMP", "Our-A2", "PhaseLift", "FontSize", 20);
ax = gca;
ax.FontSize = 20;
