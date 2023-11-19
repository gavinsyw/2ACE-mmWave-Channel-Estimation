%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% Mt: number of measurements at Tx side.
% Mr: number of measurements at Rx side.
% ULA: a structure array that groups related parameters on ULA.
% H: a structure array that groups related parameters on mmWave channel,
% please refer to Generate_Channel.m for details of its fields.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output arguments:
% F: precoder matrix (each column represents one beam pattern codebook).
% W: combiner matrix (each column represents one beam pattern codebook).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [F, W, FW] = Directional_Beam_Bayes(Mt,Mr,ULA,SNR,U,H,option)
    Nt = ULA.Nt;
    Nr = ULA.Nr;
    d = ULA.d;
    lambda = ULA.lambda;
    AoD_range = -90:1:90;
    M = Mt * Mr;
    
    if option == 1
        % generate candidate set by directional beam
        [F_try, W_try] = Directional_Beam_Angular(90,90,ULA,H);
        FW_try = kron(transpose(F_try), W_try');
    elseif option == 2
        % generate candidate set by random beam
        % generate several random beams, then select beams from them
        candidate_size = 90;
        Np = ULA.Phase_Bit^2;
        randt = randi([0 Np-1],Nt,candidate_size);
        F_try = exp(1j*randt*pi/Np)/sqrt(Nt);
        randr = randi([0 Np-1],Nr,candidate_size);
        W_try = exp(1j*randr*pi/Np)/sqrt(Nr);
        FW_try = kron(transpose(F_try), W_try');
    end
            
    noise_mtx = zeros(U, Nr*Nt, Nr*Nt);
    for i = 1:U
        for j = 1:Nr*Nt
            noise_mtx(i, j, j) = H.vecH(i,j)^-1;
        end
    end

    [fwlist, ~] = bayesAopt_complex(FW_try, M, 'K', db2pow(SNR)*squeeze(noise_mtx));

    F = zeros(Nt, Mt);
    W = zeros(U, Nr, Mr);
    FW_tmp = FW_try(fwlist, :);
    disp(["Matrix rank", rank(FW_tmp)]);
    FW = zeros(U, Mt*Mr, Nt*Nr);
    FW(1, :, :) = FW_tmp;
end