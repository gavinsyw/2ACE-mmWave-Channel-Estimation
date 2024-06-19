%% FUNCTION DESCRIPTION
% Test the ADMM method and Bayes experiment design with real RSS trace
% measured at different phase states (random phase states). After gaining
% the CSI by the proposed algorithm and baseline algorithms described in
% the MobiHoc paper, the test is hold for prediction RSS at other
% beamforming ways and evaluate their errors at predicting RSS values.

%% Initialization 
clc; clear all; close all
format long;
profile on;

%% Add paths
% In this section, three 3rd party software components are included
% 1. SparsePR: Matlab Software for Sparse Phase Retrieval
% 2. Generalized Approximate Message Passing (GAMP) in MRI MATLAB package
% 3. Orthogonal Matching Pursuit (OMP) and Compressive Sampling Matched
% Pursuit (CoSaMP) should be installed
folder_name = ...
[
"/Numerical_Simulation/3rd_software_component";
"/Numerical_Simulation/src";
"/3rd_software_component";
"/src";
];
for i=1:length(folder_name)
    Lib_path = char(strcat(pwd,folder_name(i)));
    addpath(genpath(Lib_path));
end

%% ULA parameter
ULA.lambda = 3*10^8/(60.48*10^9);  % Wavelength of 60.48 GHz signal
ULA.d = 3.055*10^(-3);             % Spacing distance between two neighboring antenna elements
ULA.Phase_Bit = 3;                 % Number of bits of phase shifters
ULA.Nt = 12;                       % Number of transmitter antenna
ULA.Nr = 12;                       % Number of receiver antenna
ULA.NQt = 4*ULA.Nt;           % Initial value of number of AoD quantization
ULA.NQr = 4*ULA.Nr;           % Initial value of number of AoA quantization

%% Setting of searching range in degree
Searching_Area = 90; 

%% Other parameters
Beampattern_Mode = 'Random_Phase_State'; % Available options shown in Generate_Sensing_Matrix.m
Add_Noise_Flag = 1;         % Whether the gaussian noise is introduced in the Rx side
on_grid = 0;                % Whether the generated AoD and AoA is on the quantized grid
L = 3;                      % Number of the dominant path
Rician_K = 5;               % Number of the other paths
Fix_angle_flag = 0;         % Whether fix the AoD and AoA for algorithm debug and test 
Show_Beam_Pattern_Flag = 0; % Whether show the exploited beam patterns in algorithm
Show_leakeage_Flag = 0;     % Whether show the quantization of the generated channel
plot_flag = 0;              % Whether show the recovery performance of the sparse channel vector
U = 1;                      % Number of Users
         
%% The different methods (algorithms) to test
Method.PhaseLift = 0;       % High computation and memory requirement
Method.CPRL = 0;            % High computation and memory requirement
Method.PRGAMP = 0;          % Bad performance by directly using PRGAMP algorithm
Method.MySparsePL = 0;      % Bad performance by directly using alternative minimization
Method.PLOMP = 0;           % Proposed algorithm by using OMP
Method.PLGAMP = 0;          % Proposed algorithm by using GAMP
Method.PerfectPhaseCS = 0;  % Perfect phase benchmark
Method.NoisyPhaseCS = 0;    % Noisy phase benchmark
Method.ADMM = 1;
Method.ADMMLowRankV1 = 1;
Method.ADMMLowRankV2 = 0;
Method.ADMMLowRankV4 = 1;
Method.Number = Method.PhaseLift + ...
                Method.CPRL + ...
                Method.PRGAMP + ...
                Method.MySparsePL + ...
                Method.PLOMP + ...
                Method.PLGAMP + ...
                Method.PerfectPhaseCS + ...
                Method.NoisyPhaseCS + ...
                Method.ADMM + ...
                Method.ADMMLowRankV1 + ...
                Method.ADMMLowRankV2 + ...
                Method.ADMMLowRankV4;
Method.State = [Method.PhaseLift ...
                Method.CPRL ...
                Method.PRGAMP ...   
                Method.MySparsePL ...
                Method.PLOMP ...
                Method.PLGAMP ...
                Method.PerfectPhaseCS ... 
                Method.NoisyPhaseCS ...
                Method.ADMM ...
                Method.ADMMLowRankV1 ...
                Method.ADMMLowRankV2 ...
                Method.ADMMLowRankV4];
Method.Name_of_method = ["PhaseLift","CPRL","PRGAMP","MySparsePL","PLOMP",...
                         "PLGAMP","PerfectPhaseCS","NoisyPhaseCS","ADMM",...
                         "ADMMLowRankV1","ADMMLowRankV2","ADMMLowRankV4"]; 
Method.Color_Line = ["gd-", "g+-", "yx-", "yo-", "r^-",...
                     "r*-","bo-","mx-","gd-", "g+-", "yx-", "yo-"];
                
%% Load data
% load rss data
load(["../Dataset/rss_trace_movement_simu_12x12_brownian/codebook_simu_movement.mat"]);

% load codebook
load(["../Dataset/rss_trace_movement_simu_12x12_brownian/rss_simu_movement.mat"]);

cb = codebook;
rss_final = db2pow(rss_trace);

% total number of time windows
Tw = size(cb, 1) / 100;

% initial number of measurement
M_max = 400;
threshold = 0.2;
for method = Method.Name_of_method
    Evaluation.(method) = [];
    M_t.(method) = zeros(Tw, 1);
end

%% Generate a sparse channel representation for arbitrary H
% Generate channel with number of path and AoA/AoD setting 
H = Generate_Channel(ULA, L, Searching_Area, Rician_K, on_grid, Fix_angle_flag);

% Generate the sparse formulation of channel
Sparse_Channel_Representation = Sparse_Channel_Formulation(ULA,L,H,Show_leakeage_Flag);
z = Sparse_Channel_Representation.z;

%% Perform Channel Estimation
n = length(z);
s = L;

% time slots are divided into 496 pieces, with 100 probes in each piece
for method = Method.Name_of_method
    if Method.(method)
        M = 80;
        beams_sliding_window = [];
        for t = 1:Tw
            disp(["t: ", t, " M: ", M]);
            M_t.(method)(t) = M;
            current_beams = 100*(t-1)+1:1:100*(t-1)+M;
            beams_sliding_window = [beams_sliding_window current_beams];
            beams_sliding_window = beams_sliding_window(max(numel(beams_sliding_window)-400, 1):numel(beams_sliding_window));

            cb_test = cb(100*(t-1)+M+1:100*t,:);
            rss_test = rss_final(100*(t-1)+M+1:100*t);

            % Generate Measurement
            measurements = rss_final(beams_sliding_window);
            beams = cb(beams_sliding_window, :);

            % Recovery
            recovered_Channel = Recover_Channel(measurements, beams, Method, ...
                    ULA, s, Sparse_Channel_Representation.AD);

            % Evaluation
            H = recovered_Channel.(method);
            [evaluation, rss_eval] = Evaluate_rss(H, cb_test, rss_test);
            Evaluation.(method)=[Evaluation.(method) evaluation];
            if evaluation < threshold
                M = max(0, M-floor(M/5)-1);
            else
                M = min(80, M+floor(M/5)+1);
            end
        end
    end
end

%% Plot result
figure();
plot(pow2db(Evaluation.ADMM(1:200)), "LineWidth", 2);
hold on;
plot(pow2db(Evaluation.ADMMLowRankV1(1:200)), "LineWidth", 2);
hold on;
plot(pow2db(Evaluation.ADMMLowRankV4(1:200)), "LineWidth", 2);
% hold on;
% plot(pow2db(Evaluation.PLOMP), 'g-^', "LineWidth", 2);
% hold on;
% plot(pow2db(Evaluation.PLGAMP), 'bx-', "LineWidth", 2);
legend("ADMM", "ADMMLowRankV1","ADMMLowRankV4", "FontSize", 12, "FontName", "Times");
xlabel("T", "FontSize", 18, "FontName", "Times");
ylabel("Error (dB)", "FontSize", 18, "FontName", "Times");
grid on;

figure();
plot(M_t.ADMM(1:200), "LineWidth", 2);
hold on;
plot(M_t.ADMMLowRankV1(1:200), "LineWidth", 2);
hold on;
plot(M_t.ADMMLowRankV4(1:200), "LineWidth", 2);
legend("ADMM", "ADMMLowRankV1", "ADMMLowRankV4", "FontSize", 12, "FontName", "Times");
xlabel("T", "FontSize", 18, "FontName", "Times");
ylabel("Probes / 100 packets", "FontSize", 18, "FontName", "Times");
grid on;
