%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2019 Yi Zhang and The University of Texas at Austin 
%  
% Permission is hereby granted, free of charge, to any person obtaining a 
% copy of this software and associated documentation files (the
% "Software"), to deal in the Software without restriction, including 
% without limitation the rights to use, copy, modify, merge, publish, 
% distribute, sublicense, and/or sell copies of the Software, and to permit
% persons to whom the Software is furnished to do so, subject to the 
% following conditions:
% 
% The above copyright notice and this permission notice shall be included
% in all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
% NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
% DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
% OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
% USE OR OTHER DEALINGS IN THE SOFTWARE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If you use this code or any (modified) part of it in any publication,
% please cite:
%
% Yi Zhang, Kartik Patel, Sanjay Shakkottai, and Robert W. Heath Jr.. 2019. 
% Side-information-aided Non-coherent Beam Alignment Design for Millimeter 
% Wave Systems. In MobiHoc '19: The Twentieth ACM International Symposium 
% on Mobile Ad Hoc Networking and Computing, July 02-05, 2019, Catania, 
% Italy. ACM, New York, NY, USA, 10 pages.
%
% Author: Yi Zhang
% Contact email: yi.zhang.cn@utexas.edu 
% Last modified: Apr. 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function description:
% This function manages all the recovery (compressive phase retrieval) 
% algorithms and the related benchmarking algorithms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% measurements: a structure array that groups the different types of  
% received noised measurements. Please refer to Generate_Measurement.m for 
% the details of its fields.
% measurementMat: sensing matrix, it is stored in the field.
% "measurementMat" of the output of the function Generate_Sensing_Matrix
% z: vectorized and quantized mmWave channel. it is stored in the field "z" 
% of the output of the function Sparse_Channel_Formulation.
% s: level of sparsity, i.e. number of large dominant paths. 
% plot_flag: whether show the recovery performance of the sparse channel 
% vector.
% noise_power: power of noise.
% Method: a structure array that groups related parameters on the recovery
% algorithms that needs to be tested.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output arguments:
% recoveredSig_Set: a structure array that groups all the recovered sparse
% signal (estimated z) by different methods. Each field of recoveredSig_Set
% is a column vector.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function recoveredSig_Set = MyCPR_H(measurements, Sensing_Matrix, z, s, plot_flag, noise_power, Method, ULA)
%% Read parameters
measurementMat = Sensing_Matrix.measurementMat;
F = Sensing_Matrix.F;
FW = Sensing_Matrix.FW;

%% Problem parameters
measurements_norm_square = measurements.measurements_norm_square;
measurements_with_perfect_phase = measurements.measurements_with_perfect_phase;
measurements_with_noisy_phase = measurements.measurements_with_noisy_phase;
TrueSig = z;
[U, m, n] = size(measurementMat);
if plot_flag % Print out the key parameters of the problem to be solved
    fprintf('Problem dimension - %d\n', n);
    fprintf('Sparsity - %d\n', s );
    fprintf('No. of measurements_norm_square - %d\n', m);
end


%% PhaseLift
if Method.PhaseLift
    recoveredSig_Set.PhaseLift = zeros(n,U);
    for i = 1:U
        PhaseLift_singleRx = MyPhaseLift(measurements_norm_square(:,i), squeeze(measurementMat(i,:,:)));
        if plot_flag
            Plot_Recovery_Performance(TrueSig(i,:),recoveredSig_Set.PhaseLiftsingleRx,'PhaseLift');
        end
        recoveredSig_Set.PhaseLift(:,i) = PhaseLift_singleRx;
    end
end


%% CPRL
if Method.CPRL
    recoveredSig_Set.CPRL = zeros(n,U);
    for i = 1:U
        CPRL_singleRx = MyCPRL(measurements_norm_square(:,i), squeeze(measurementMat(i,:,:)));
        if plot_flag
            Plot_Recovery_Performance(TrueSig(i,:),CPRL_singleRx,'CPRL');
        end
        recoveredSig_Set.CPRL(:,i) = CPRL_singleRx;
    end
end


% ----------------------------------------------------------------
% ---------------- NOT AVAILABLE CURRENTLY -----------------------
% ----------------------------------------------------------------
% %% PR-GAMP
% if Method.PRGAMP
%     opt_pr.xreal = 0;
%     opt_pr.xnonneg = 0;
%     opt_pr.tuneOn = 1;
%     opt_pr.sparseRat = s/n;
%     opt_pr.SNRdBinit = 20;
%     opt_pr.SNRdB = 25;
%     opt_pr.verbose = 0;
%     opt_pr.maxTry = 100;
%     opt_pr.nitEM = 50;     
%     recoveredSig_Set.PRGAMP = MyPRGAMP(measurements_norm_square, measurementMat, opt_pr);
%     if plot_flag
%         Plot_Recovery_Performance(TrueSig,recoveredSig_Set.PRGAMP,'PRGAMP');
%     end
% end


% %% SparsePhaseLift
% if Method.MySparsePL
%     opt_AltMin.MaxIte = 1e4;
%     opt_AltMin.s = s;
%     opt_AltMin.opt_epslion = 0.001;
%     recoveredSig_Set.MySparsePL = ...
%         MySparsePL(measurements_norm_square, measurementMat, opt_AltMin);
%     if plot_flag
%         Plot_Recovery_Performance(TrueSig,recoveredSig_Set.MySparsePL,'MySparsePL');
%     end
% end


%% Two_Stage_Recovery (Our proposed algorithm)
if Method.PLOMP || Method.PLGAMP
    recoveredSig_Set.PLOMP = zeros(n, U);
    recoveredSig_Set.PLGAMP = zeros(n, U);
    for i = 1:U
        [plomp_singleRx, plgamp_singleRx] = ...
            My_TwoStage_Recovery(measurements_norm_square(:,i), ...
            squeeze(measurementMat(i,:,:)), s, noise_power, plot_flag, Method);
        recoveredSig_Set.PLOMP(:, i) = plomp_singleRx;
        recoveredSig_Set.PLGAMP(:, i) = plgamp_singleRx;
    end
%     if plot_flag && Method.PLOMP
%         Plot_Recovery_Performance(TrueSig,recoveredSig_Set.PLOMP,'PLOMP');
%     end
%     if plot_flag && Method.PLGAMP
%         Plot_Recovery_Performance(TrueSig,recoveredSig_Set.PLGAMP,'PLGAMP');
%     end
end


%% Conventional CS via measurement with perfect phase information
if Method.PerfectPhaseCS
    [w, ~] = size(F);
    recoveredSig_Set.PerfectPhaseCS = zeros(n,U);
    for i = 1:U
        PerfectPhaseCS_singleRx = My_Conventional_CS(measurements_with_perfect_phase(:,i), squeeze(measurementMat(i,:,:)), s, noise_power);
%         PerfectPhaseCS_singleRx = My_Unconventional_CS(transpose(measurements_with_perfect_phase(i,:)), F);
        if plot_flag
            Plot_Recovery_Performance(TrueSig(i,:),PerfectPhaseCS_singleRx,'PerfectPhaseCS');
        end
        recoveredSig_Set.PerfectPhaseCS(:,i) = PerfectPhaseCS_singleRx;
    end
end


%% Conventional CS via measurement with noisy phase information
if Method.NoisyPhaseCS
    [w, ~] = size(F);
    recoveredSig_Set.NoisyPhaseCS = zeros(n, U);
    for i = 1:U
        NoisyPhaseCS_singleRx = My_Conventional_CS(measurements_with_noisy_phase(:,i), squeeze(measurementMat(i,:,:)), s, noise_power);
%         NoisyPhaseCS_singleRx = My_Unconventional_CS(transpose(measurements_with_noisy_phase(i,:)), F);
        if plot_flag
            Plot_Recovery_Performance(TrueSig(i,:),recoveredSig_Set.NoisyPhaseCS,'NoiseyPhaseCS');
        end
        recoveredSig_Set.NoisyPhaseCS(:,i) = NoisyPhaseCS_singleRx;
    end
end

%% ADMM channel matrix recovery with signal recovery
if Method.ADMM
    [~, ~, w] = size(FW);
    recoveredSig_Set.ADMM = zeros(w, U);
    for i = 1:U
        version = 0;
        [ADMM_singleRx, C, D] = ADMM_v2(sqrt(measurements_norm_square(:,i)), squeeze(FW(i,:,:)), ULA.Nt, ULA.Nr, version);
        recoveredSig_Set.ADMM(:, i) = ADMM_singleRx;
    end
end

%% ADMM method with low rank assumption
if Method.ADMMLowRank
    [~, ~, w] = size(FW);
    recoveredSig_Set.ADMMLowRank = zeros(w, U);
    for i = 1:U
        version = 1;
        [ADMM_singleRx, C, D] = ADMM_v2(sqrt(measurements_norm_square(:,i)), squeeze(FW(i,:,:)), ULA.Nt, ULA.Nr, version);
        recoveredSig_Set.ADMMLowRank(:, i) = ADMM_singleRx;
    end
end

%% ADMM method with low rank assumption
if Method.ADMMLowRankV2
    [~, ~, w] = size(FW);
    recoveredSig_Set.ADMMLowRankV2 = zeros(w, U);
    for i = 1:U
        version = 2;
        [ADMM_singleRx, C, D] = ADMM_v2(sqrt(measurements_norm_square(:,i)), squeeze(FW(i,:,:)), ULA.Nt, ULA.Nr, version);
        recoveredSig_Set.ADMMLowRankV2(:, i) = ADMM_singleRx;
    end
end