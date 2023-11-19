function [final_group, calibrated_phase_shift, calibration_factor, calibration_factor_digit] = group_ant_kernel(phaseshift, TX, group_size, phase_bit, azimuth_aoa, elevation_aoa, freq)

% greedy heuristic for antenna grouping.
% use the first available antenna as the reference antenna, and then select
% the antennas that has the least phase difference between it as the group
% members.

phaseshift_vec =  phaseshift(TX);
available_phase = (0:1/(2^phase_bit):1) * 2 * pi;

locality_phase_azimuth = cos(azimuth_aoa) * cos(elevation_aoa) * 2 * 0.58 * pi;
locality_phase_elevation = cos(elevation_aoa) * sin(azimuth_aoa) * 2 * 0.58 * pi;

group_num = length(TX) / group_size;
final_group = zeros(group_num, group_size);
allocated_flag = zeros(1, length(TX));
calibrated_phase_shift = zeros(1, 32);
calibration_factor = zeros(1, 32);

for i = 1:group_num
    available_elements = find(allocated_flag==0);
    first_element = min(available_elements);
    cost_vec = zeros(1, numel(available_elements)-1);
    selected_phaseshift = zeros(1, numel(available_elements)-1);
    [first_azimuth_phaseoffset, first_elevation_phaseoffset] = get_location_phaseoffset(locality_phase_azimuth, locality_phase_elevation, TX(first_element));
    for j = 2:numel(available_elements)
        current_element = available_elements(j);
        [cur_azimuth_phaseoffset, cur_elevation_phaseoffset] = get_location_phaseoffset(locality_phase_azimuth, locality_phase_elevation, TX(current_element));
        azimuth_phaseoffset = cur_azimuth_phaseoffset - first_azimuth_phaseoffset;
        elevation_phaseoffset = cur_elevation_phaseoffset - first_elevation_phaseoffset;
        total_phase_offset = mod(phaseshift_vec(current_element)-phaseshift_vec(first_element)+azimuth_phaseoffset+elevation_phaseoffset, 2*pi);
%         total_phase_offset = cur_total_phaseoffset - first_total_phaseoffset;
%         total_phase_offset = mod(phaseshift_vec(current_element)-phaseshift_vec(first_element)+total_phase_offset, 2*pi);
        [cost_vec(j-1), selected_phaseshift(j-1)] = min(abs(total_phase_offset-available_phase));
    end
    [~, idx] = sort(cost_vec, "ascend");
    selected_elements = [first_element, available_elements(1+idx(1:group_size-1))];
    selected_phaseshift = selected_phaseshift(idx(1:group_size-1));
    
    allocated_flag(selected_elements) = 1;
    
    final_group(i, :) = selected_elements;
    for j = 2:group_size
%         current_element = selected_elements(j);
%         [cur_azimuth_phaseoffset, cur_elevation_phaseoffset] = get_location_phaseoffset(locality_phase_azimuth, locality_phase_elevation, current_element);
%         azimuth_phaseoffset = cur_azimuth_phaseoffset - first_azimuth_phaseoffset;
%         elevation_phaseoffset = cur_elevation_phaseoffset - first_elevation_phaseoffset;
%         total_phase_offset = mod(phaseshift_vec(first_element)-phaseshift_vec(current_element)+azimuth_phaseoffset+elevation_phaseoffset, 2*pi);
%         [~, idx] = min([abs(available_phase(selected_phaseshift(j-1))-total_phase_offset), mod(abs(available_phase(selected_phaseshift(j-1))+total_phase_offset), 2*pi)]);
%         if idx == 1
        calibrated_phase_shift(TX(selected_elements(j))) = mod(phaseshift_vec(selected_elements(j)) - available_phase(selected_phaseshift(j-1)), 2*pi);
        calibration_factor(TX(selected_elements(j))) = available_phase(selected_phaseshift(j-1));
        calibration_factor_digit = round(calibration_factor./(pi/2));
        calibration_factor_digit(calibration_factor_digit ==4) =0;
%         else
%             calibrated_phase_shift(selected_elements(j)) = mod(phaseshift_vec(selected_elements(j))) - available_phase(selected_phaseshift(j-1));
%         end
    end
end
final_group = TX(final_group);
end

function [azimuth_phaseoffset, elevation_phaseoffset] = get_location_phaseoffset(locality_phase_azimuth, locality_phase_elevation, element)
%      ground_truth_x_idx = [0,0,-1,-1,1,1,1,0,1,-1,0,-1,0,1,0,1,3,3,4,4,2,2,2,3,4,2,3,4,3,2,3,2];
%      ground_truth_y_idx = [0,1,0,1,1,-1,0,-1,2,3,2,2,3,4,4,3,0,1,0,1,1,-1,0,-1,3,2,2,2,3,4,4,3];  
    %ground_truth_x_idx = [0,0,1,1,-1,-1,-1,0,-1,1,0,1,0,-1,0,-1,-3,-3,-4,-4,-2,-2,-2,-3,-4,-2,-3,-4,-3,-2,-3,-2];
    %ground_truth_y_idx = [0,-1,0,-1,-1,1,0,1,-2,-3,-2,-2,-3,-4,-4,-3,0,-1,0,-1,-1,1,0,1,-3,-2,-2,-2,-3,-4,-4,-3];
    ground_truth_x_idx = [0,0,1,1,1,1,1,0,1,1,0,1,0,1,0,1,3,3,4,4,2,2,2,3,4,2,3,4,3,2,3,2];
    ground_truth_y_idx = [0,1,0,1,1,1,0,1,2,3,2,2,3,4,4,3,0,1,0,1,1,1,0,1,3,2,2,2,3,4,4,3];
    x_idx = ground_truth_x_idx(element);
    y_idx = ground_truth_y_idx(element);
    azimuth_phaseoffset = locality_phase_azimuth * x_idx;
    elevation_phaseoffset = locality_phase_elevation * y_idx;
end

