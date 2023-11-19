function H2 = Generate_Dynamic_Channel(H0, max_angle_change, ULA, L, Rician_K, U)

lambda = ULA.lambda;
d = ULA.d;
Nt = ULA.Nt;
Nr = ULA.Nr; 

H2.AoD_array = zeros(U, L);
H2.AoA_array = zeros(U, L);
for i = 1:U
    H2.AoD_array(i, :) = H0.AoD_array(i,:) + rand(1, L)*max_angle_change;
    H2.AoA_array(i, :) = H0.AoA_array(i,:) + rand(1, L)*max_angle_change;
end

H2.AoD_array_radian = deg2rad(H2.AoD_array);
H2.AoA_array_radian = deg2rad(H2.AoA_array);
H2.h_array = zeros(U, L);
for i = 1:U
%     H2.h_array(i, :) = (randn(1,L)+randn(1,L)*1j)*1/sqrt(2); 
    H2.h_array(i, :) = H0.h_array(i, :);
    H2.h_array = H2.h_array/norm(H2.h_array); % Normalized complex channel gain
end
% Rician channel model
if L > 1
    H2.Rician_K = 0;
else
    H2.Rician_K = Rician_K;
end
H2.NLOS_h_array = zeros(U, H2.Rician_K);
AoD_Rician_Radian = zeros(U, H2.Rician_K);
AoA_Rician_Radian = zeros(U, H2.Rician_K);
for i = 1:U
    H2.NLOS_h_array(i, :) = (randn(1,H2.Rician_K)+randn(1,H2.Rician_K)*1j)*1/sqrt(2);
    H2.NLOS_h_array(i, :) = H2.NLOS_h_array/norm(H2.NLOS_h_array);
    AoD_Rician_Radian(i, :) = unifrnd(-pi/2,pi/2,1,H2.Rician_K);
    AoA_Rician_Radian(i, :) = unifrnd(-pi/2,pi/2,1,H2.Rician_K);
end

% Dominant Channel Matrix
H2.H_Matrix_Dominant = zeros(U, Nr, Nt);
for i = 1:U
    ATx = zeros(Nt,L);
    ARx = zeros(Nr,L);
    for l=1:1:L
        ATx(:,l) = 1/sqrt(Nt)*transpose(exp(-1i*2*pi/lambda*d*sin(H2.AoD_array_radian(i, l)).*(0:1:Nt-1)));
        ARx(:,l) = 1/sqrt(Nr)*transpose(exp(-1i*2*pi/lambda*d*sin(H2.AoA_array_radian(i, l)).*(0:1:Nr-1)));
    end
    H2.H_Matrix_Dominant(i, :, :) = sqrt(Nt*Nr)*ARx*diag(H2.h_array(i, :))*ATx';    
end

% Undominant Channel Matrix
H2.H_Matrix_Undominant = zeros(U, Nr, Nt);
for i = 1:U
    ATx = zeros(Nt,H2.Rician_K);
    ARx = zeros(Nr,H2.Rician_K);
    for k=1:1:H2.Rician_K
        ATx(:,k) = 1/sqrt(Nt)*transpose(exp(-1i*2*pi/lambda*d*sin(AoD_Rician_Radian(i, k)).*(0:1:Nt-1)));
        ARx(:,k) = 1/sqrt(Nr)*transpose(exp(-1i*2*pi/lambda*d*sin(AoA_Rician_Radian(i, k)).*(0:1:Nr-1)));
    end
    H2.H_Matrix_Undominant(i, :, :) = sqrt(Nt*Nr)*ARx*diag(H2.NLOS_h_array(i, :))*ATx';    
end

% Channel Matrix
K_Factor_dB = 7;
K_Factor = 10^(K_Factor_dB/10);
if H2.Rician_K ~= 0
    H2.H_Matrix = sqrt(K_Factor/(K_Factor+1))*H2.H_Matrix_Dominant + sqrt(1/(K_Factor+1))*H2.H_Matrix_Undominant;
else
    H2.H_Matrix = H2.H_Matrix_Dominant;
end
H2.vecH = zeros(U, Nt*Nr);
for i = 1:U
    H2.vecH(i, :) = vec(H2.H_Matrix(i, :));  
end
H2.AoD_range = H0.AoD_range;
H2.AoA_range = H0.AoA_range;


end