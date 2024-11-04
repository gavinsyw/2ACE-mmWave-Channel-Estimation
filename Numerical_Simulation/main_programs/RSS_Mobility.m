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
ULA.Nt = 16;                       % Number of transmitter antenna
ULA.Nr = 16;                       % Number of receiver antenna
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
load(["../Dataset/random_probe_5m_indoor/random_probe_16ant_codebook.mat"]);

% load rss
load(["../Dataset/random_probe_rss_LOS_indoor_complex_moving_040421.mat"]);

% ant_idx = [21,5,26,9]; % 4 ant
% ant_idx = [23,7,21,5,26,9,32,16]; % 8 ant
ant_idx = [17,23,7,1,18,21,5,2,27,26,9,11,29,32,16,13];   % 16 ant
% ant_idx = 1:32;
cb = kron(rx_codebook(:,ant_idx), tx_codebook(:,ant_idx));

% cb = codebook;
% rss_final = rss_trace;

rss_fct = 1e6;

% total number of time windows
Tw = 100;
T_size = 62;

cb = cb(1:1:T_size*Tw, :);
% cb = repmat(cb, Tw, 1);
rss_final = db2pow(rss(1:Tw, 1:T_size)') * rss_fct;
rss_final = vec(rss_final);

% initial number of measurement
M_max = 50;
Mw_max = 80;
threshold = 0.3;
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
        H = zeros(ULA.Nt*ULA.Nr, 1);
        M = 0;
        beams_sliding_window = [];
        for t = 1:Tw
            disp(["t: ", t, " M: ", M]);
            M_t.(method)(t) = M;
            current_beams = T_size*(t-1)+1:1:T_size*t;

            cb_current = cb(T_size*(t-1)+1:T_size*t,:);
            rss_current = rss_final(T_size*(t-1)+1:T_size*t);
            
            % evaluate the performance using the previous inferred channel
            % using the current beams
            [evaluation, rss_eval] = Evaluate_rss(H, cb_current, rss_current);
            Evaluation.(method)=[Evaluation.(method) evaluation];
            if evaluation < threshold 
                M = 0;
                beams_sliding_window = [beams_sliding_window current_beams];
                beams_sliding_window = beams_sliding_window(max(length(beams_sliding_window)-Mw_max, 1) : length(beams_sliding_window));
            else
                M = min(ceil(M*1.2+1), Mw_max);
                beams_sliding_window = [beams_sliding_window current_beams];
                beams_sliding_window = beams_sliding_window(max(length(beams_sliding_window)-Mw_max, 1) : length(beams_sliding_window));
            end
            
            % Generate Measurement
            measurements = rss_final(beams_sliding_window);
            beams = cb(beams_sliding_window, :);

            % Recovery
            recovered_Channel = Recover_Channel(measurements, beams, Method, ...
                    ULA, s, Sparse_Channel_Representation.AD);

            % Evaluation
            H = recovered_Channel.(method);

        end
    end
end

%% save result
save("mobility_real_trace.mat");

%% read csv
csv_data = csvread("../Dataset/Untitled from 2021-04-04 18_59_39 +0000.csv",1,4, [1 4 5719 4]);
% plot(csv_data);
% xlabel("t", 'Interpreter', "latex");
% ylabel("X value", "Interpreter", "latex");

%% Plot result
plot_init;
load("mobility_real_trace.mat");
figure();

hl1 = line(1:Tw, pow2db(Evaluation.ADMMLowRankV4), 'Color','r');
xlabel("T", "Interpreter", "latex");
ylabel("Error (dB)", "Interpreter", "latex");
ax1 = gca;
set(ax1,'XColor','r','YColor','r');

ax2 = axes('Position',get(ax1,'Position'),...
           'XAxisLocation','top',...
           'YAxisLocation','right',...
           'Color','none',...
           'XColor','k','YColor','k');
hl2 = line(1:length(csv_data), csv_data, 'Color','k','Parent',ax2);
axis([0, length(csv_data), -inf, inf]);
xlabel("t", "Interpreter", "latex");
ylabel("X", "Interpreter", "latex");
grid on;

figure();
% yyaxis left;
hl3 = line(1:Tw, M_t.ADMMLowRankV4/T_size * 100, "LineWidth", 2);
xlabel("T", "Interpreter", "latex");
ylabel("Probes / 100 packets", "Interpreter", "latex");
ax1 = gca;
set(ax1,'XColor','r','YColor','r');

ax2 = axes('Position',get(ax1,'Position'),...
           'XAxisLocation','top',...
           'YAxisLocation','right',...
           'Color','none',...
           'XColor','k','YColor','k');
hl4 = line(1:length(csv_data), csv_data, 'Color','k','Parent',ax2);
axis([0, length(csv_data), -inf, inf]);
xlabel("t", "Interpreter", "latex");
ylabel("X", "Interpreter", "Latex");
grid on;
