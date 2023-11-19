function [H_amp,H_angle] = channel_recovery_ADMM_v2_simulation_directional(tx_ant_num, rx_ant_num, cb_amp, cb_angle, rss_final, seed_id)
    %% FUNCTION DESCRIPTION
    % Test the ADMM method and Bayes experiment design with real RSS trace
    % measured at different phase states (random phase states). After gaining
    % the CSI by the proposed algorithm and baseline algorithms described in
    % the MobiHoc paper, the test is hold for prediction RSS at other
    % beamforming ways and evaluate their errors at predicting RSS values.

    %% Initialization 
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
    "/3rd_software_component";
    "/src";
    ];
    for i=1:length(folder_name)
        Lib_path = char(strcat(pwd,folder_name(i)));
        addpath(genpath(Lib_path));
    end

    %% ULA parameter
    tx_ant_num = double(tx_ant_num);
    rx_ant_num = double(rx_ant_num);
    ULA.lambda = 3*10^8/(60.48*10^9);  % Wavelength of 60.48 GHz signal
    ULA.d = 2.9*10^(-3);             % Spacing distance between two neighboring antenna elements
    ULA.Phase_Bit = 2;                 % Number of bits of phase shifters
    ULA.Nt = double(tx_ant_num);                       % Number of transmitter antenna
    ULA.Nr = double(rx_ant_num);                       % Number of receiver antenna
    ULA.NQt = 4*ULA.Nt;           % Initial value of number of AoD quantization
    ULA.NQr = 4*ULA.Nr;           % Initial value of number of AoA quantization

    % number of measurement
    %M = double(M);

    %% Setting of searching range in degree
    Searching_Area = 180; 

    %% Other parameters
    Beampattern_Mode = 'Directional_Beam_Angular'; % Available options shown in Generate_Sensing_Matrix.m
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
    Method.PRGAMP = 0;          % Bad performance by directly us           ing PRGAMP algorithm
    Method.MySparsePL = 0;      % Bad performance by directly using alternative minimization
    Method.PLOMP = 1;           % Proposed algorithm by using OMP
    Method.PLGAMP = 1;          % Proposed algorithm by using GAMP
    Method.ADMM = 0;            % ADMM
    Method.ADMMLowRankV1 = 0;     % Low-rank ADMM
    Method.ADMMLowRankV2 = 0;     % Low-rank ADMM
    Method.ADMMLowRankV3 = 0;     % Low-rank ADMM
    Method.ADMMLowRankV4 = 0;     % Low-rank ADMM
    Method.Number = Method.PhaseLift + ...
                    Method.CPRL + ...
                    Method.PRGAMP + ...
                    Method.MySparsePL + ...
                    Method.PLOMP + ...
                    Method.PLGAMP + ...
                    Method.ADMM + ...
                    Method.ADMMLowRankV1 + ...
                    Method.ADMMLowRankV2 + ...
                    Method.ADMMLowRankV3 + ...
                    Method.ADMMLowRankV4;
    Method.State = [Method.PhaseLift ...
                    Method.CPRL ...
                    Method.PRGAMP ...   
                    Method.MySparsePL ...
                    Method.PLOMP ...
                    Method.PLGAMP ...
                    Method.ADMM ...
                    Method.ADMMLowRankV1 ...
                    Method.ADMMLowRankV2...
                    Method.ADMMLowRankV3...
                    Method.ADMMLowRankV4];
    Method.Name_of_method = ["PhaseLift","CPRL","PRGAMP","MySparsePL","PLOMP",...
                             "PLGAMP","ADMM","ADMMLowRankV1","ADMMLowRankV2","ADMMLowRankV3", "ADMMLowRankV4"]; 
    
    rng(4096);
    
    if (tx_ant_num == 4)||(rx_ant_num == 4)
        %load('./codebook/random_probe_4ant_codebook.mat');
        %ant_idx = [5, 9, 21, 26]; % 4 ant
        %M = [8, 16, 24, 32, 40, 48, 56, 64];
    elseif (tx_ant_num == 8)||(rx_ant_num == 8)
        Mt = round(linspace(2,sqrt(4*tx_ant_num*rx_ant_num),8));
        Mr = Mt;
    elseif (tx_ant_num == 16)||(rx_ant_num == 16)
        Mt = round(linspace(2,sqrt(4*tx_ant_num*rx_ant_num),8));
        Mr = Mt;
    elseif (tx_ant_num == 17)||(rx_ant_num == 17)
        Mt = round(linspace(2,sqrt(4*(tx_ant_num-1)*(rx_ant_num-1)),8));
        Mr = Mt;
    elseif (tx_ant_num == 32)||(rx_ant_num == 32)
        Mt = round(linspace(2,sqrt(4*tx_ant_num*rx_ant_num),8));
        Mr = Mt;
    elseif (tx_ant_num == 36)||(rx_ant_num == 36)
        Mt = round(linspace(2,sqrt(4*tx_ant_num*rx_ant_num),8));
        Mr = Mt;
    else
        error('Number of antenna on Tx and Rx must be 4/8/16/32!')
    end

    cb = cb_amp.*exp(1j*cb_angle);

    %% Split Training and Testing Dataset
    % assign the first 49000 entries as training dataset. Beams are selected
    % among them. Then the last 600 entries act as the testing dataset.
    % Evaluation are done among the 600 entries to test the errors of RSS.

    % add a multiplication factor for rss to move it near 1.
    rss_fct = 1e5/3;

    H_out = zeros(length(Mt), Method.Number, tx_ant_num*rx_ant_num);
    for i = 1: length(Mt)

        M_cur = Mt(i);
        if M_cur <= 32
            indexing = round(linspace(1,32,M_cur));
        else
            indexing = round(linspace(1,32,32));
        end

        cb_train = cb(indexing,indexing,:);

        beams = reshape(cb_train,[M_cur^2, tx_ant_num*rx_ant_num]);

        rss_train = rss_final(indexing,indexing);
        rss_train = reshape(rss_train, [M_cur^2, 1]);
        measurements = sqrt(db2pow(rss_train) / 1000) * rss_fct;

        %% Generate a sparse channel representation for arbitrary H
        % Generate channel with number of path and AoA/AoD setting 
        H = Generate_Channel(ULA, L, Searching_Area, Rician_K, on_grid, Fix_angle_flag, U);
        
        % Generate the sparse formulation of channel
        Sparse_Channel_Representation = Sparse_Channel_Formulation(ULA,L,H,Show_leakeage_Flag,U);
        %z = Sparse_Channel_Representation.z;

        %% Perform Channel Estimation
        %n = length(z);
        s = L;

        % Generate Sensing Matrix 
        %picked_beams =...
        %    Generate_Sensing_Matrix(Beampattern_Mode,M_cur,ULA,cb_train,rss_train);
        % Sensing_Matrix =...
        %     Generate_Sensing_Matrix(Beampattern_Mode,Mt,Mr,ULA,Sparse_Channel_Representation,Show_Beam_Pattern_Flag,L,H)

        % Recovery
        recovered_Channel = Recover_Channel(measurements, beams, Method, ...
                ULA, s, Sparse_Channel_Representation.AD);
        
        ind = 1;
        for method = Method.Name_of_method 
            if Method.(method)
                disp(method);
                H_cur = recovered_Channel.(method);
                H_out(i,ind,:) = H_cur./rss_fct; %sum(H_cur,1);
                ind = ind + 1;
            end
        end
        disp(['finished measurement with ' num2str(M_cur^2) ' probes']);
    end
    H_out(isnan(H_out)) = 0;
    H_amp = double(abs(H_out));
    H_angle = double(angle(H_out));
end
