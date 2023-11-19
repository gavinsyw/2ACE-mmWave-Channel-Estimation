function [X, Y, converged] = ADMM_v2(measurements, FW, TX, RX, version)

%   seed = 5;
%   rng(seed);
  
  % get input
  [T, ~] = size(FW);
  B = measurements;
  A = FW;
%   rz = min([T TX RX]);
  rz = 3;
  
  % Generate input
%   A = randn(T, TX*RX) + 1j*randn(T, TX*RX);
%   rz = 2;
%   Z_gt = (randn(TX,rz) + 1j*randn(TX, rz)) * (randn(rz, RX) + 1j*randn(rz ,RX));
%   X_gt = reshape(Z_gt, TX*RX, 1);
%   Y_gt = A * X_gt;
%   B = abs(Y_gt);
  
  % Given A and c, solve x
  if version == 0
    [X, Y, converged] = inferMinL2(A, B);
  elseif version == 1
    [X, Y, converged] = inferLowRank(A, B, TX, RX);
  elseif version == 2
    [X, Y, converged] = inferLowRankV2(A, B, TX, RX);
  elseif version == 3
    [X, Y, converged] = inferLowRankV3(A, B, TX, RX);
  elseif version == 4
    %[X, Y, converged] = inferLowRankV4(A, B, TX, RX);
    [X, Y, converged] = inferLowRankV4_multi(A, B, TX, RX);
  else
    R = TX;
    RZ = rz + 2;
    for iter = 1:3
      [X, Y, converged] = inferLowRankV2(A, B, TX, RX, RZ, R);
      if converged
        break
      else
        R = R + floor(TX/2);
        RZ = RZ + 2;
      end
    end
  end
  
%   X_ratio = X./X_gt;
%   error = norm(X_ratio - mean(X_ratio)*ones(size(X_ratio)))/norm(X_ratio);
