clear
clc

load('hardware_phaseoffset.mat');

ele = 90/180*pi;
azi = 90/180*pi; %[45, 60, 75, 90,105,120,135]
freq = 60.48e9;

ant_choice = [1,2,3,4,5,6,7,8,17,18,19,20,21,22,23,24];


[final_group_4, calibrated_phase_shift_4, calibration_factor_4, calibration_factor_digit_4] = group_ant_kernel(antenna_phase_shifts, ant_choice, 4, 2, ele, azi);
[final_group_2, calibrated_phase_shift_2, calibration_factor_2, calibration_factor_digit_2] = group_ant_kernel(antenna_phase_shifts, ant_choice, 2, 2, ele, azi);

antenna_phase_shifts_bit = round(-antenna_phase_shifts./((pi)/2));
antenna_phase_shifts_bit(antenna_phase_shifts_bit<0) = antenna_phase_shifts_bit(antenna_phase_shifts_bit<0) + 4;
antenna_phase_shifts_bit(antenna_phase_shifts_bit == 4) = 0;
antenna_phase_shifts_bit = antenna_phase_shifts_bit.';
antenna_phase_shifts_bit([9,10,11,12,13,14,15,16,25,26,27,28,29,30,31,32]) = 0;
antenna_phase_shifts_bit(ant_choice)

% for i = 1:4
%     for j = 1:4
%         %calibrated_phase_shift_4(final_group_4(i,j)) = calibrated_phase_shift_4(final_group_4(i,j)) - antenna_phase_shifts(final_group_4(i,1));
%         calibration_factor_digit_4(final_group_4(i,j)) = round((calibration_factor_digit_4(final_group_4(i,j))*pi/2 - antenna_phase_shifts(final_group_4(i,1)))./(pi/2));
%         calibration_factor_digit_4(calibration_factor_digit_4 < 0) =calibration_factor_digit_4(calibration_factor_digit_4 < 0) + 4;
%                 calibration_factor_digit_4(calibration_factor_digit_4 > 4) =calibration_factor_digit_4(calibration_factor_digit_4 > 4) - 4;
%         calibration_factor_digit_4(calibration_factor_digit_4 ==4) =0;
%     end
% end
% 
% for i = 1:8
%     for j = 1:2
%         %calibrated_phase_shift_2(final_group_2(i,j)) = calibrated_phase_shift_2(final_group_2(i,j)) - antenna_phase_shifts(final_group_2(i,1));
%         calibration_factor_digit_2(final_group_2(i,j)) = round((calibration_factor_digit_2(final_group_2(i,j))*pi/2 - antenna_phase_shifts(final_group_2(i,1)))./(pi/2));
%         calibration_factor_digit_2(calibration_factor_digit_2 < 0) =calibration_factor_digit_2(calibration_factor_digit_2 < 0) + 4;
%         calibration_factor_digit_2(calibration_factor_digit_2 > 4) =calibration_factor_digit_2(calibration_factor_digit_2 > 4) - 4;
%         calibration_factor_digit_2(calibration_factor_digit_2 ==4) =0;
%     end
% end