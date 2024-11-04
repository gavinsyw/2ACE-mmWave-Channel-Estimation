function Sensing_Matrix = Generate_random_beam(Mt,Mr,ULA,...
            Sparse_Channel_Representation,Show_Beam_Pattern_Flag,L,H,seed)

%% Parameter fetching
lambda = ULA.lambda;
d = ULA.d;
Nt = ULA.Nt;
Nr = ULA.Nr;
Phase_Bit = ULA.Phase_Bit;        
AD = Sparse_Channel_Representation.AD;

% generate random numbers for each measurement one by one
% this ensures that for as Mt increases, the set of measurements
% contains all measurements for smaller Mt
rng(seed);
Np = Phase_Bit^2;
randt = zeros(Nt,Mt);
RandRs = zeros(Nr,Mr);

for k = 1:Mt
  randt(:,k) = randi([0 Np-1],Nt,1);
  RandRs(:,k) = randi([0 Np-1],Nr,1);
end

% sensing matrix at transmitter side
F = exp(1j*randt*2*pi/Np)/sqrt(Nt);

% sensing matrix at receiver side
W = F;

for k = 1:Mt*Mr
    randrt(:,k) = randi([0 Np-1],Nt*Nr,1);
end
FW = transpose(exp(1j*randrt*2*pi/Np))/sqrt(Nt*Nr);

A = FW*AD;

Sensing_Matrix.F = F;
Sensing_Matrix.W = W;  
Sensing_Matrix.FW = FW;
Sensing_Matrix.AD = AD;        
Sensing_Matrix.measurementMat = A;    
end