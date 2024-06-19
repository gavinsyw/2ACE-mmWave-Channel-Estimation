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
ULA.Phase_Bit = 2;                 % Number of bits of phase shifters
ULA.Nt = 32;                       % Number of transmitter antenna
ULA.Nr = 32;                       % Number of receiver antenna
ULA.NQt = 4*ULA.Nt;           % Initial value of number of AoD quantization
ULA.NQr = 4*ULA.Nr;           % Initial value of number of AoA quantization

%% Number of measurement range 
% M = [ 5  10  15  20  25];               % Number of measurements at Tx and Rx side
M = 200:1:200;
% G = [25 40 55 60 70];               % Number of AoD and AoA quantization

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
         
%% The different methods (algorithms) to test
Method.PhaseLift = 0;       % High computation and memory requirement
Method.CPRL = 0;            % High computation and memory requirement
Method.PRGAMP = 0;          % Bad performance by directly using PRGAMP algorithm
Method.MySparsePL = 0;      % Bad performance by directly using alternative minimization
Method.PLOMP = 0;           % Proposed algorithm by using OMP
Method.PLGAMP = 0;          % Proposed algorithm by using GAMP
Method.PerfectPhaseCS = 0;  % Perfect phase benchmark
Method.NoisyPhaseCS = 0;    % Noisy phase benchmark
Method.ADMM = 0;
Method.ADMMLowRankV1 = 0;
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
                 
for method = Method.Name_of_method
    Evaluation.(method) = [];
end

%% Load data
% load codebook
load(["../Dataset/random_probe_5m_indoor/random_probe_32ant_codebook.mat"]);

% load rss data
rss_data_path = "random_probe_rss_5m_32ant_LOS_reshaped.mat";
load(["../Dataset/random_probe_5m_indoor/"+rss_data_path]);

% ant_idx = [21,5,26,9]; % 4 ant
% ant_idx = [23,7,21,5,26,9,32,16]; % 8 ant
% ant_idx = [17,23,7,1,18,21,5,2,27,26,9,11,29,32,16,13];   % 16 ant
ant_idx = 1:32;

cb = kron(rx_codebook(:,ant_idx), tx_codebook(:,ant_idx));
rss_final = rss_out;

% cb = codebook;
% rss_final = rss_trace;

%% Split Training and Testing Dataset
% assign the first 49000 entries as training dataset. Beams are selected
% among them. Then the last 600 entries act as the testing dataset.
% Evaluation are done among the 600 entries to test the errors of RSS.

% add a multiplication factor for rss to move it near 1.
rss_fct = 1e7;

record_H = zeros(30, 32, 32);

%% Generate a sparse channel representation for arbitrary H
% Generate channel with number of path and AoA/AoD setting 
H = Generate_Channel(ULA, L, Searching_Area, Rician_K, on_grid, Fix_angle_flag);

% Generate the sparse formulation of channel
Sparse_Channel_Representation = Sparse_Channel_Formulation(ULA,L,H,Show_leakeage_Flag);
z = Sparse_Channel_Representation.z;

%% Perform Channel Estimation
n = length(z);
s = L;

for i = 1:30
    cb_train = cb((i-1)*200+1:i*200, :);
    rss_train = db2pow(rss_final((i-1)*200+1:i*200, 1)) * rss_fct;


    % Generate Sensing Matrix 
    picked_beams = 1:200;
    beams = cb_train(picked_beams, :);

    % Generate Measurement
    measurements = rss_train(picked_beams);

    % Recovery
    recovered_Channel = Recover_Channel(measurements, beams, Method, ...
            ULA, s, Sparse_Channel_Representation.AD);

    for method = Method.Name_of_method
        if Method.(method)
            % Evaluation
            H_inferred = recovered_Channel.ADMMLowRankV4;
            record_H(i,:,:) = reshape(H_inferred, 32, 32);
        end
    end

end

%% save result
save("Inferred_channel_"+rss_data_path, "record_H");