function ideal_steering_vector = IdealSteeringVectorPerPanel(azim_range, elev_range, antenna_map, beam_map)
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
    
    antenna_number = length(antenna_map);
    antenna_geometry = zeros(3, antenna_number);
    antenna_spacing = 0.58;
    for ii = 1:antenna_number
        y_id = floor((ii-1)/6);
        x_id = mod((ii-1),6);
        antenna_geometry(:,ii) = [x_id; y_id; 0] * antenna_spacing;
    end
    
    ideal_steering_vector = zeros(antenna_number, azim_num, elev_num);
    for ii = 1:azim_num
        for jj = 1:elev_num
            ideal_steering_vector(:,ii,jj) = exp(1j * 2 * pi * aoa_unit(:,ii,jj).' * antenna_geometry)';
        end
    end
    
    ideal_steering_vector = ideal_steering_vector(beam_map,:,:);
    for ii = 32:-1:1
        ideal_steering_vector(ii,:,:) = ideal_steering_vector(ii,:,:) .* conj(ideal_steering_vector(1,:,:));
    end
end