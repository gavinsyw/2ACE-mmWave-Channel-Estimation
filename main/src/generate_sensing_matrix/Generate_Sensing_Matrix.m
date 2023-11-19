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
% This function generates the sensing matrix.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% Method: a structure array that groups related parameters on the
% activation setting of the different methods (algorithms) to test.
% Mt: number of measurements at Tx side.
% Mr: number of measurements at Rx side.
% ULA: a structure array that groups related parameters on ULA.
% Sparse_Channel_Representation: a structure array that groups related 
% parameters on quantized sparse mmWave channel. Please refer to 
% Sparse_Channel_Formulation.m for details of its fields.
% Show_Beam_Pattern_Flag: whether to plot the beam pattern associated with
% the generated sensing matrix.
% L: Number of the dominant path.
% H: a structure array that groups related parameters on mmWave channel,
% please refer to Generate_Channel.m for details of its fields.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output arguments:
% Sensing_Matrix: a structure array that groups related parameters sensing
% matrix, please refer to Generate_Sensing_Matrix.m for details of its
% fields.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Sensing_Matrix = Generate_Sensing_Matrix(Method,Mt,Mr,ULA,Sparse_Channel_Representation,Show_Beam_Pattern_Flag,L,H,U,SNR,seed)
    %% Parameter fetching
    lambda = ULA.lambda;
    d = ULA.d;
    Nt = ULA.Nt;
    Nr = ULA.Nr;
    Phase_Bit = ULA.Phase_Bit;
    AD = Sparse_Channel_Representation.AD;
    [~, AD_len] = size(AD);
    
    %% Generate beam patterns 
    switch Method
        case 'Random_Phase_State'                   % random gaussian sensing matrix without constraints on phase array 
            % generate random numbers for each measurement one by one
            % this ensures that for as Mt increases, the set of measurements
            % contains all measurements for smaller Mt
            rng(seed);
            Np = Phase_Bit^2;
            randt = zeros(Nt,Mt);
            RandRs = zeros(U,Nr,Mr);

            for k = 1:Mt
              randt(:,k) = randi([0 Np-1],Nt,1);
              for i = 1:U
                RandRs(i,:,k) = randi([0 Np-1],Nr,1);
              end
            end

            % sensing matrix at transmitter side
            F = exp(1j*randt*2*pi/Np)/sqrt(Nt);
            
            % sensing matrix at receiver side
            W = zeros(U, Nr, Mr);
            FW = zeros(U, Mt*Mr, Nt*Nr);
            A = zeros(U, Mt*Mr, AD_len);
            
            for i = 1:U
                randrt = zeros(Nt*Nr,Mt*Mr);
                for k = 1:Mt*Mr
                    randrt(:,k) = randi([0 Np-1],Nt*Nr,1);
                end
                randr = squeeze(RandRs(i,:,:));
                W_single = exp(1j*randr*2*pi/Np)/sqrt(Nr);
%                 W_single = 1;
%                 W(i,:,:) = W_single;
                FW(i,:,:) = transpose(exp(1j*randrt*2*pi/Np))/sqrt(Nt*Nr);
%                 FW(i,:,:) = kron(transpose(F),W_single');
                A(i,:,:) = squeeze(FW(i,:,:))*AD;
            end

%         ----------------------------------------------------------------
%         ---------------- NOT AVAILABLE CURRENTLY -----------------------
%         ----------------------------------------------------------------
%         case 'Gaussian_Infinite_Bits_Phase'         % random gaussian sensing matrix with infinite-bit phase shifters   
%             F = (randn(Nt,Mt) + 1i*randn(Nt,Mt));
%             F = exp(1j*angle(F))/sqrt(Nt);
%             W = (randn(Nr,Mr) + 1i*randn(Nr,Mr));
%             W = exp(1j*angle(W))/sqrt(Nr);
% 
%         case 'Gaussian_Finite_Bits_Phase'           % random gaussian sensing matrix with finite-bit phase shifters
%             F = (randn(Nt,Mt) + 1i*randn(Nt,Mt));
%             F = Quantize_PS(F,Phase_Bit);
%             W = (randn(Nr,Mr) + 1i*randn(Nr,Mr));
%             W = Quantize_PS(W,Phase_Bit);
% 
%         case 'Gaussian_Finite_Bits_Phase_Low_Rank'  % random gaussian sensing matrix with finite-bit phase shifters and with rank reduced operation
%             F = (randn(Nt,Mt) + 1i*randn(Nt,Mt));
%             F = Quantize_PS(F,Phase_Bit);
%             rk = round(Mt*0.4);
%             for i=1:1:rk
%                 coef = randn(Mt-rk,1);
%                 coef = coef/norm(coef);
%                 F(:,i) = F(:,rk+1:Mt)*coef;
%             end
%             %F = Quantize_PS(F,Phase_Bit);
%             W = (randn(Nr,Mr) + 1i*randn(Nr,Mr));
%             W = Quantize_PS(W,Phase_Bit);    
%             rk = round(Mr*0.4);
%             for i=1:1:rk
%                 coef = randn(Mt-rk,1);
%                 coef = coef/norm(coef);
%                 W(:,i) = W(:,rk+1:Mr)*coef;
%             end            
%             %W = Quantize_PS(W,Phase_Bit);
%             FW = kron(transpose(F),W');
%             [ui,~]=size(unique(FW,'row'));
%             fprintf('Rank of FW %d \n\n', rank(FW));
%             fprintf('Unique row of FW %d \n\n', ui);
%        
%         case 'Region_Random_Beam'                   % random beam pattern within a certain region
%             [F, W] = Region_Random_Beam(Mt,Mr,ULA,H);
%             
%         case 'Directional_Random_Beam'              % directional beam pattern with random gain in spatial domain
%             [F, W] = Directional_Random_Beam(Mt,Mr,ULA,H);
%         ----------------------------------------------------------------

        case 'Directional_Beam'                     % directional beam pattern with uniform gain in spatial domain
            Rank_Eliminated = 0;
            [F, W_single] = Directional_Beam(Mt,Mr,ULA,H,Rank_Eliminated);
            W = zeros(U, Nr, Mr);
            FW = zeros(U, Mt*Mr, Nt*Nr);
            A = zeros(U, Mt*Mr, AD_len);
            for i = 1:U
                W(i,:,:) = W_single;
                FW(i,:,:) = kron(transpose(F),W_single');
                A(i,:,:) = squeeze(FW(i,:,:))*AD;
            end
            
        case 'Directional_Beam_Angular'             % directional beam pattern with uniform gain in angular domain
            [F, W_single] = Directional_Beam_Angular(Mt,Mr,ULA,H);   
            W = zeros(U, Nr, Mr);
            FW = zeros(U, Mt*Mr, Nt*Nr);
            A = zeros(U, Mt*Mr, AD_len);
            for i = 1:U
                W(i,:,:) = W_single;
                FW(i,:,:) = kron(transpose(F),W_single');
                A(i,:,:) = squeeze(FW(i,:,:))*AD;
            end
            
%         ----------------------------------------------------------------
%         ---------------- NOT AVAILABLE CURRENTLY -----------------------
%         ----------------------------------------------------------------
%         case 'Directional_Beam_Low_Rank'            % directional beam pattern with uniform gain in spatial domain and with reduced rank operation
%             Rank_Eliminated = 2;
%             [F, W] = Directional_Beam(Mt,Mr,ULA,H,Rank_Eliminated);
%                
%         case 'ProtoType_Pattern'                    % using real beam pattern that is used by the phased array hardware
%             load('prmCode.mat');
%             load('BeamformerQuant_F.mat');          % This mat file stores the beamformer
%             F = BeamformerQuant_F(:,:,prmCode.Ind);
%             load('BeamformerQuant_W.mat');          % This mat file stores the combiner
%             W = BeamformerQuant_W(:,:,prmCode.Ind);
%         ----------------------------------------------------------------

        case 'Directional_Beam_Bayes'
            option = 1;
            A = zeros(U, Mt*Mr, AD_len);
            [F, W, FW] = Directional_Beam_Bayes(Mt,Mr,ULA,SNR,U,H,option);
%             for i = 1:U
%                 A(i,:,:) = squeeze(FW(i,:,:))*AD;
%             end
        
        case 'Random_Beam_Bayes'
            option = 2;
            A = zeros(U, Mt*Mr, AD_len);
            [F, W, FW] = Directional_Beam_Bayes(Mt,Mr,ULA,SNR,U,H,option);
        
        case 'Directional_Beam_Bayes_v2'
            W = zeros(U, Nr, Mr);
            FW = zeros(U, Mt*Mr, Nt*Nr);
            A = zeros(U, Mt*Mr, AD_len);
            [F, W_single] = Directional_Beam_Bayes_v2(Mt,Mr,ULA,AD,SNR,U,H);
            for i = 1:U
                W(i,:,:) = W_single;
                FW(i,:,:) = kron(transpose(F), W_single');
                A(i,:,:) = squeeze(FW(i,:,:))*AD;
            end
        otherwise
            ;
    end
    if Show_Beam_Pattern_Flag
        show_beam_pattern(lambda, d, F, Method);
%         show_beam_pattern(lambda, d, W);
    end

%     ----------------------------------------------------------------
%     ---------------- NOT AVAILABLE CURRENTLY -----------------------
%     ----------------------------------------------------------------
%     %% Output result
%     if rank(F) < min(Mt,Nt)
%         fprintf('Repeat beam pattern with rank of F is %d \n', rank(F));
%     end
%     if rank(W) < min(Mr,Nr)
%         fprintf('Repeat beam pattern with rank of W is %d \n', rank(W));
%     end        

    Sensing_Matrix.F = F;
    Sensing_Matrix.W = W;  
    Sensing_Matrix.FW = FW;
    Sensing_Matrix.AD = AD;        
    Sensing_Matrix.measurementMat = A;    
end


