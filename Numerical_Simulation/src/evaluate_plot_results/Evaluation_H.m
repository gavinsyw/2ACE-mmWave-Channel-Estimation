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
% This function helps to plot the simulation results of the proposed
% algorithm and the related benchmarking algorithms.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% recoveredSig: X-axis label (character array).
% H: a structure array that groups related parameters on mmWave channel,
% please refer to Generate_Channel.m for details of its fields.
% Sparse_Channel_Representation: a structure array that groups related 
% parameters on quantized sparse mmWave channel. Please refer to
% Sparse_Channel_Formulation.m for the details of its fields.
% ULA: a structure array that groups the related parameters on ULA.
% SNR: SNR in dB.
% Mtr: a structure array that groups Mt and Mr, where Mt denotes the
% number.
% of measurements used in the Tx side while Mr denotes the number of
% measurements used in the Rx side.
% noise_power: power of noise.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output arguments:
% Evaluation_result: a structure array that groups the results with 
% different evaluation metrics.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [Evaluation_result, H_Estimated] = Evaluation_H(recoveredSig,L,H,Sparse_Channel_Representation,ULA,SNR,Mtr,noise_power)
    % load parameters
    Nt = ULA.Nt;
    Nr = ULA.Nr;
    Phase_Bit = ULA.Phase_Bit;
    
    % estimate normalize error
    vecH_Estimated = recoveredSig;
    phaseFac = exp( 1i* angle( (vecH_Estimated'*H.vecH)/(H.vecH'*H.vecH) ) );
    vecH_Estimated = vecH_Estimated.*phaseFac;
    H_estimated = reshape(vecH_Estimated,Nr,Nt);
    H_real = H.H_Matrix;
    MSE_H = norm(H_estimated-H_real,'fro')^2/norm(H_real,'fro')^2;
    
    X_gt = H.vecH;
    X = vecH_Estimated;
    MSE_H = norm(X_gt - (X'*X_gt)/(X'*X)*X,'fro')^2/norm(X_gt,'fro')^2;
    
%     disp(MSE_H);
    
    % beamforming gain
    H_Estimated = reshape(vecH_Estimated, Nr, Nt);
    [U,~,V] = svd(H_Estimated);
    w_Unconstrained = U(:,1);
    f_Unconstrained = V(:,1);
    w_PerCSI = Quantize_PS(w_Unconstrained,Phase_Bit);
    f_PerCSI = Quantize_PS(f_Unconstrained,Phase_Bit);
    w_dig = w_Unconstrained ./ norm(w_Unconstrained);
    f_dig = f_Unconstrained ./ norm(f_Unconstrained);
    
    gain_ana = abs(w_PerCSI' * squeeze(H.H_Matrix) * f_PerCSI);
    gain_dig = abs(w_dig' * squeeze(H.H_Matrix) * f_dig);
    
    % 1D projection error
    [U,S,V] = svd(H.H_Matrix);
    X_gt = U(:,1)*S(1,1)*V(:,1)';
    X_gt = X_gt(:);

    [U,S,V] = svd(H_Estimated);
    X = U(:,1)*S(1,1)*V(:,1)';
    X = X(:);

    proj_error = norm(X_gt - (X'*X_gt)/(X'*X)*X)/norm(X_gt);

    
    Evaluation_result=[0;...
                   0;...
                   0;...
                   MSE_H;
                   0;
                   0;
                   0;
                   0;
                   0;
                   0;
                   0;
                   0;
                   gain_ana;
                   gain_dig;
                   proj_error];
%     Evaluation_result=[Num_Quantization_Error;...
%                        AoDA_Err;...
%                        MSE_ARM;...
%                        MSE_H;...
%                        Rate_Unconstrained;...
%                        Rate_PerCSI;...
%                        Rate_Est;...
%                        Rate_Est_AoDA;...
%                        Rate_BSweep1;...
%                        Rate_BSweep2;...
%                        AoDA_Err_Sweep1;...
%                        AoDA_Err_Sweep2];
%      save AoDA_Err AoDA_Err;
end

