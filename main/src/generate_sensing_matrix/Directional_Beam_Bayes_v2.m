%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% Mt: number of measurements at Tx side.
% Mr: number of measurements at Rx side.
% ULA: a structure array that groups related parameters on ULA.
% H: a structure array that groups related parameters on mmWave channel,
% SNR: global SNR readings from probing signals,
% U: total number of receivers,
% please refer to Generate_Channel.m for details of its fields.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output arguments:
% F: precoder matrix (each column represents one beam pattern codebook).
% W: combiner matrix (each column represents one beam pattern codebook).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function usage:
% Generate sensing matrix at the transmitter and receiver side by
% Bayes-A optimization. The transmitter antennas are equivalently separated
% into U parts, and then adapt bayes-opt function to find the values for
% beamforming.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [F, W] = Directional_Beam_Bayes_v2(Mt,Mr,ULA,AD,SNR,U,H)
    Nt = ULA.Nt;
    Nr = ULA.Nr;
    d = ULA.d;
    lambda = ULA.lambda;
%     AoD_range = H.AoD_range(1):1:H.AoD_range(2);
    AoD_range = -90:1:90;
    max_ind = length(AoD_range);
    F_try = zeros(Nt,max_ind);
    W_try = zeros(Nr,max_ind);
    n_sep = floor(Nt / U);  % number of antennas in each group
    for j = 1:U
        for i=1:max_ind
            F_try(n_sep*(j-1)+1:n_sep*j,i) = transpose(exp(-1i*2*pi*d/lambda*sind(AoD_range(i)).*(0:1:n_sep-1)))/sqrt(Nt);
        end
    end
    for i = 1:max_ind
        W_try(:,i) = transpose(exp(-1i*2*pi*d/lambda*sind(AoD_range(i)).*(0:1:Nr-1)))/sqrt(Nr);
    end

    M = Mt * Mr;
%     FW_try = kron(transpose(F_try), W_try')*AD;
    [~,n] = size(AD);
    F_try_new = F_try * sqrt(db2pow(SNR));
    K_all = zeros(U, n, n);
    for i = 1:U
        H_mtx = reshape(H.H_Matrix(i,:,:), Nr, Nt);
        K_all(i, :, :) = find_K(H_mtx, n, ULA);
    end
    [fwlist, ~] = MyBayesAopt(transpose(F_try_new)*AD, M, U, 'K', K_all);
%     [flist, ~] = bayesAopt(F_try', Mt, 'K', abs(H.H_Matrix).^-1);
%     [wlist, ~] = bayesAopt(W_try', Mr, 'K', abs(H.H_Matrix).^-1);
%     flist = floor(fwlist./max_ind+1);
%     wlist = mod(fwlist, max_ind)+1;
    
    % take the most commonly appeared numbers
%     % flist
%     b = sort(flist);
%     [~, indx] = unique(b);
%     appear_num = [indx(2:length(indx)); length(b)+1] - indx;
%     [~, y] = sort(appear_num, "descend");
%     flist = flist(y(1:Mt));
%     % same for wlist
%     b = sort(wlist);
%     [~, indx] = unique(b);
%     appear_num = [indx(2:length(indx)); length(b)+1] - indx;
%     [~, y] = sort(appear_num, "descend");
%     wlist = wlist(y(1:Mr));

    F_try_quant = Quantize_PS(F_try,ULA.Phase_Bit);
%     W_try_quant = Quantize_PS(W_try,ULA.Phase_Bit);
    F = F_try_quant(:, fwlist);
%     W = W_try_quant(:, wlist);
    W = 1;
end