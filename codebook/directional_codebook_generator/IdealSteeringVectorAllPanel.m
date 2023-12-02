%#############################################################################
%Project 2ACE
%PI: Prof. Lili Qiu @ UT Austin/MSRA
%Participants: Yiwen Song @ CMU, Changhan Ge @ UT Austin, Yin Zhang
%Copyright @ The University of Texas at Austin, 2023
%Partially inherited from Kun Qian@University of Virginia & Xinyu Zhang @ UCSD
%#############################################################################

function ideal_steering_vector = IdealSteeringVectorAllPanel(azim_range, elev_range, antenna_map, beam_map, panel_map, module_map)
    ideal_steering_vector_per_panel = IdealSteeringVectorPerPanel(azim_range, elev_range, antenna_map, beam_map);
    
    azim_num = length(azim_range);
    elev_num = length(elev_range);
    aoa_unit = zeros(3, azim_num, elev_num);
    for ii = 1:azim_num
        for jj = 1:elev_num
            aoa_unit(:,ii,jj) = [cosd(azim_range(ii))*cosd(elev_range(jj)); ...
                sind(azim_range(ii))*cosd(elev_range(jj)); ...
                sind(elev_range(jj))];
        end
    end
    
    beam_number = length(beam_map);
    panel_number = length(panel_map);
    panel_geometry = zeros(3, panel_number);
    panel_spacing = 0.58*6;
    for ii = 1:panel_number
        y_id = floor((ii-1)/2);
        x_id = mod((ii-1),2);
        panel_geometry(:,ii) = [x_id; y_id; 0] * panel_spacing;
    end
    
    ideal_steering_vector_panel = zeros(panel_number, azim_num, elev_num);
    for ii = 1:azim_num
        for jj = 1:elev_num
            ideal_steering_vector_panel(:,ii,jj) = exp(1j * 2 * pi * aoa_unit(:,ii,jj).' * panel_geometry)';
        end
    end
    
    ideal_steering_vector_panel = ideal_steering_vector_panel(module_map,:,:);
    
    ideal_steering_vector = zeros(panel_number, beam_number, azim_num, elev_num);
    for ii = 1:panel_number
        ideal_steering_vector(ii,:,:,:) = ideal_steering_vector_per_panel .* repmat(ideal_steering_vector_panel(ii,:,:), beam_number, 1, 1);
    end
end