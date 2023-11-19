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
% This function generates the sparse mmWave channel.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% ULA: a structure array that groups related parameters on ULA.
% L: number of the dominant path.
% Searching_Area: searching range in degree.
% Rician_K: number of the other paths.
% on_grid: whether the generated AoD and AoA is on the quantized grid.
% Fix_angle_flag: whether fixing AoD and AoA for algorithm debug and test.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output arguments:
% H: a structure array that groups related parameters on mmWave channel.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function H = Generate_Channel(ULA, L, Searching_Area, Rician_K, on_grid, Fix_angle_flag, U)
    lambda = ULA.lambda;
    d = ULA.d;
    Nt = ULA.Nt;
    Nr = ULA.Nr; 
    NQt = ULA.NQt;
    NQr = ULA.NQr; 
    AoD_range = [-Searching_Area/2 Searching_Area/2];
    AoA_range = AoD_range; 
    if Fix_angle_flag == 1
        H.AoD_array = zeros(U, 1);
        H.AoA_array = ones(U, 1) * 15;
    else
        H.AoD_array = zeros(U, L);
        H.AoA_array = zeros(U, L);
        for i = 1:U
            pretect_zone = 0;
            H.AoD_array(i, :) = unifrnd(AoD_range(1)+pretect_zone,AoD_range(2)-pretect_zone,1,L);
            H.AoA_array(i, :) = unifrnd(AoA_range(1)+pretect_zone,AoA_range(2)-pretect_zone,1,L);
        end
    end
    if on_grid
        % Quantization of Spatial Domain
        partition_T = linspace(-1,1,NQt+1);
        partition_T = partition_T(1:end-1);
        partition_R = linspace(-1,1,NQr+1);
        partition_R = partition_R(1:end-1);
        AoD_Virtual_quantized_array = partition_T;
        AoA_Virtual_quantized_array = partition_R;
        for i = 1:U
            for l=1:1:L
                [~, p] = min(abs(AoD_Virtual_quantized_array - sind(H.AoD_array(i, l))));
                H.AoD_array(i, l) = asind(AoD_Virtual_quantized_array(p));
                [~, p] = min(abs(AoA_Virtual_quantized_array - sind(H.AoA_array(i, l))));
                H.AoA_array(i, l) = asind(AoA_Virtual_quantized_array(p));
            end
        end
    end
    H.AoD_array_radian = deg2rad(H.AoD_array);
    H.AoA_array_radian = deg2rad(H.AoA_array);
    H.h_array = zeros(U, L);
    for i = 1:U
        H.h_array(i, :) = (randn(1,L)+randn(1,L)*1j)*1/sqrt(2); 
        H.h_array = H.h_array/norm(H.h_array); % Normalized complex channel gain
    end
    % Rician channel model
    if L > 1
        H.Rician_K = 0;
    else
        H.Rician_K = Rician_K;
    end
    H.NLOS_h_array = zeros(U, H.Rician_K);
    AoD_Rician_Radian = zeros(U, H.Rician_K);
    AoA_Rician_Radian = zeros(U, H.Rician_K);
    for i = 1:U
        H.NLOS_h_array(i, :) = (randn(1,H.Rician_K)+randn(1,H.Rician_K)*1j)*1/sqrt(2);
        H.NLOS_h_array(i, :) = H.NLOS_h_array/norm(H.NLOS_h_array);
        %H.NLOS_h_array = H.NLOS_h_array/norm(H.NLOS_h_array)*sqrt(Nt*Nr); % Normalized complex channel gain
        AoD_Rician_Radian(i, :) = unifrnd(-pi/2,pi/2,1,H.Rician_K);
        AoA_Rician_Radian(i, :) = unifrnd(-pi/2,pi/2,1,H.Rician_K);
    end

    % Dominant Channel Matrix
    H.H_Matrix_Dominant = zeros(U, Nr, Nt);
    for i = 1:U
        ATx = zeros(Nt,L);
        ARx = zeros(Nr,L);
        for l=1:1:L
            ATx(:,l) = 1/sqrt(Nt)*transpose(exp(-1i*2*pi/lambda*d*sin(H.AoD_array_radian(i, l)).*(0:1:Nt-1)));
            ARx(:,l) = 1/sqrt(Nr)*transpose(exp(-1i*2*pi/lambda*d*sin(H.AoA_array_radian(i, l)).*(0:1:Nr-1)));
        end
        H.H_Matrix_Dominant(i, :, :) = sqrt(Nt*Nr)*ARx*diag(H.h_array(i, :))*ATx';    
    end
    
    % Undominant Channel Matrix
    H.H_Matrix_Undominant = zeros(U, Nr, Nt);
    for i = 1:U
        ATx = zeros(Nt,H.Rician_K);
        ARx = zeros(Nr,H.Rician_K);
        for k=1:1:H.Rician_K
            ATx(:,k) = 1/sqrt(Nt)*transpose(exp(-1i*2*pi/lambda*d*sin(AoD_Rician_Radian(i, k)).*(0:1:Nt-1)));
            ARx(:,k) = 1/sqrt(Nr)*transpose(exp(-1i*2*pi/lambda*d*sin(AoA_Rician_Radian(i, k)).*(0:1:Nr-1)));
        end
        H.H_Matrix_Undominant(i, :, :) = sqrt(Nt*Nr)*ARx*diag(H.NLOS_h_array(i, :))*ATx';    
    end
       
    % Channel Matrix
    K_Factor_dB = 7;
    K_Factor = 10^(K_Factor_dB/10);
    if H.Rician_K ~= 0
        H.H_Matrix = sqrt(K_Factor/(K_Factor+1))*H.H_Matrix_Dominant + sqrt(1/(K_Factor+1))*H.H_Matrix_Undominant;
    else
        H.H_Matrix = H.H_Matrix_Dominant;
    end
    H.vecH = zeros(U, Nt*Nr);
    for i = 1:U
        H.vecH(i, :) = vec(H.H_Matrix(i, :));  
    end
    H.AoD_range = AoD_range;
    H.AoA_range = AoA_range;
end

