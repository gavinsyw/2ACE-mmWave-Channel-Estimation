%###########################################################################
%Project 2ACE
%PI: Prof. Lili Qiu @ UT Austin/MSRA
%Participants: Yiwen Song @ CMU, Changhan Ge @ UT Austin, Yin Zhang
%Copyright @ The University of Texas at Austin, 2023
%Partially inherited from Kun Qian@UVa & Xinyu Zhang @ UCSD
%###########################################################################

function antenna_phase_shifts = AntennaPhaseShifts(steering_phase, azim_range, elev_range, antenna_map, beam_map)
    real_steering_vector = exp(1j * steering_phase);
    ideal_steering_vector = IdealSteeringVectorPerPanel(azim_range, elev_range, antenna_map, beam_map);
    steering_vector_diff = real_steering_vector .* conj(ideal_steering_vector);
    antenna_phase_shifts = exp(1j * angle(sum(sum(steering_vector_diff,3),2)));
end