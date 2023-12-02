%###########################################################################
%Project 2ACE
%PI: Prof. Lili Qiu @ UT Austin/MSRA
%Participants: Yiwen Song @ CMU, Changhan Ge @ UT Austin, Yin Zhang
%Copyright @ The University of Texas at Austin, 2023
%Partially inherited from Kun Qian@UVa & Xinyu Zhang @ UCSD
%###########################################################################

clear;clc;
close all;

load('steering_vector_calib.mat');
antenna_phase_shifts = AntennaPhaseShifts(steering_phase, azim_range, elev_range, antenna_map, beam_map);

azim_range = [0:36:180];
elev_range = [0:36:180];
ideal_steering_vector = IdealSteeringVectorAllPanel(azim_range, elev_range, antenna_map, beam_map, panel_map, module_map);
real_steering_vector = ideal_steering_vector .* permute(repmat(antenna_phase_shifts, 1, length(panel_map), length(azim_range), length(elev_range)), [2,1,3,4]);
real_steering_phase = mod(floor((angle(real_steering_vector) + pi/4) / (pi/2)), 4);

codebook = squeeze(real_steering_phase(1,:,:,:));

[num_ant, num_cb] = size(codebook);

fprintf('[');
for i = 1:num_cb
    fprintf('''')
    fprintf('%d', codebook(:,i));
    fprintf('''')
    if i< num_cb
        fprintf(',\n');
    end
end
fprintf(']');
