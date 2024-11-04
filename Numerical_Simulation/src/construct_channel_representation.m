function H = construct_channel_representation(real_trace, L, ULA)

d = ULA.d;
Nt = ULA.Nt;
Nr = ULA.Nr; 

% still under construction
% ---------------------------
H.AoD_array = zeros(1,L);
H.AoA_array = zeros(1,L);

H.AoD_array_radian = deg2rad(H.AoD_array);
H.AoA_array_radian = deg2rad(H.AoA_array);
H.h_array = zeros(1,L);
% ---------------------------

H.H_Matrix = real_trace / norm(squeeze(real_trace), "fro");
% H.H_Matrix = real_trace;
H.H_Matrix_Dominant = H.H_Matrix;
H.H_Matrix_Undominant = zeros(size(H.H_Matrix));

H.vecH = vec(H.H_Matrix);  

% FoV AoD & AoA search range
H.AoD_range = [-90 90];
H.AoA_range = [-90 90];

end