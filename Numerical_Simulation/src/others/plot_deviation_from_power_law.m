% synthetic data
TX = 36;
RX = 36;
% H = zeros(5, TX, RX);
% for rz = 2:6
% H(rz-1,:,:) = (randn(TX,rz) + 1j*randn(TX, rz)) * (randn(rz, RX) + 1j*randn(rz, RX));
% end

% Wireless Insite data
H = load("../Dataset/H_livingroom_180Rx.mat");
H = H.H;
% H = H(1:5, :, :);

plot_flag = 0;
err = variance_with_K_singular_values(H, plot_flag);

% Our profile
f0 = 0.8;
f1 = 0.9;
f2 = 0.95;
f3 = 0.995;
std_err_A1 = zeros(1, TX);
std_err_A1(ceil(1.5*sqrt(TX)):TX) = f2;
std_err_A2 = zeros(1, TX);
% std_err_A2(ceil(0.5*sqrt(TX)):TX) = f0;
std_err_A2(ceil(0.7*sqrt(TX)):TX) = f0;
std_err_A2(ceil(sqrt(TX)):TX) = f1;
std_err_A2(ceil(1.5*sqrt(TX)):TX) = f2;
std_err_A2(ceil(3*sqrt(TX)):TX) = f3;

% deviation from our profile
deviation_A1 = zeros(1, 180);
deviation_A2 = zeros(1, 180);
for i = 1:180
    s = std_err_A1 - (1-err(i, :));
    deviation_A1(i) = sum(s(s>0)) / norm(std_err_A1);
    s = std_err_A2 - (1-err(i, :));
    deviation_A2(i) = sum(s(s>0)) / norm(std_err_A2)/3;
end

figure();
cdfplot(deviation_A1);
hold on;
cdfplot(deviation_A2);
legend("A1", "A2");



