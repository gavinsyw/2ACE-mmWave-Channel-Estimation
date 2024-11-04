% synthetic
% H = zeros(5, TX, RX);
% for rz = 2:6
% H(rz-1,:,:) = (randn(TX,rz) + 1j*randn(TX, rz)) * (randn(rz, RX) + 1j*randn(rz, RX));
% end

% wireless insite
H = load("../Dataset/H_livingroom_180Rx.mat");
H = H.H;
H = H(1:5, :, :);

l1norm = zeros(5, 1);

for i = 1:5
    H_sub = squeeze(H(i,:,:));
%     H_sub = H_sub / norm(H_sub, "fro");
    l1norm(i) = norm(H_sub(:), 1);
end

plot(l1norm);