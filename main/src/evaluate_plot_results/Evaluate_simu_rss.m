function Evaluation_result = Evaluate_simu_rss(recoveredSig,L,H,Sparse_Channel_Representation,ULA,SNR,Mtr,noise_power,U)

%% Recover H first
% load parameters
Nt = ULA.Nt;
Nr = ULA.Nr;
Phase_Bit = ULA.Phase_Bit;

vecH_Estimated = recoveredSig;
phaseFac = exp( 1i* angle( (vecH_Estimated'*transpose(H.vecH(1,:)))/(H.vecH(1,:)*H.vecH(1,:)') ) );
vecH_Estimated = vecH_Estimated*phaseFac;

%% Generate some random Beams
Mt = 10;
Mr = 10;
Np = Phase_Bit^2;
randt = randi([0 Np-1],Nt,Mt);
F = exp(1j*randt*pi/Np)/sqrt(Nt);
randr = randi([0 Np-1],Nr,Mr);
W_single = exp(1j*randr*pi/Np)/sqrt(Nr);
FW = kron(transpose(F),W_single');

%% Use the random beam to test difference
rss_estimated = FW * vecH_Estimated;
rss_gt = FW * transpose(H.vecH);
mse = mean(abs(rss_estimated - rss_gt));
mse = norm(H.vecH-vecH_Estimated) / norm(H.vecH);
H_estimated = reshape(vecH_Estimated, Nr, Nt);
H_real = reshape(H.H_Matrix, Nr, Nt);
mse = norm(H_estimated-H_real,'fro')^2/norm(H_real,'fro')^2;

%% Calculate signal strength according to optimal analog/digital beamforming
H_Estimated = reshape(vecH_Estimated, Nr, Nt);
[U2,~,V] = svd(H_Estimated);
w_Unconstrained = exp(1j*angle(U2(:,1)))./sqrt(Nr);
f_Unconstrained = exp(1j*angle(V(:,1)))./sqrt(Nt);
w_PerCSI = Quantize_PS(w_Unconstrained,Phase_Bit);
f_PerCSI = Quantize_PS(f_Unconstrained,Phase_Bit);
sig_strength_ana = abs(w_PerCSI' * squeeze(H.H_Matrix) * f_PerCSI);
sig_strength_dig = abs(w_Unconstrained' * squeeze(H.H_Matrix) * f_Unconstrained);

%% combine the results
Evaluation_result = [mse; sig_strength_ana; sig_strength_dig];
end