function [X, Y, quality] = inferLowRankV3(A, B, tx, rx, lambda, r, mu0, rho, cc_frac, tol_rel, tol_abs, maxiter)
  if nargin < 5, lambda = 0; end
  if nargin < 6, r = 20; end
  if nargin < 7, mu0 = 0.001; end
  if nargin < 8, rho = 1.03; end
  if nargin < 9, cc_frac = 0.95; end 
  if nargin < 10, tol_rel = 1e-4; end
  if nargin < 11, tol_abs = 1e-8; end
  if nargin < 12, maxiter = 500; end
  
  [m,n] = size(A);
  r = min([r m n]);
  
  %----------------------------------------------------------
  % pre-processing
  %----------------------------------------------------------
  % scale A so norm(A,'fro') = sqrt(m)
  % scale B so norm(B) = 1
  
  A_norm = norm(A,'fro') / sqrt(m);
  if (A_norm < tol_abs)
    A_norm = 1;
  end

  B_norm = norm(B);
  if (B_norm < tol_abs), 
    B_norm = 1; 
  end
  
  A = A / A_norm;
  B = B / B_norm;

  %----------------------------------------------------------
  % solve the same problem using 95% constraints and check for consistency on the rest 
  %----------------------------------------------------------
  [m,n] = size(A);
  train_idx = randsample(m, floor(m*cc_frac));
  test_idx = setdiff(1:m, train_idx);
  A_train = A(train_idx,:);
  B_train = B(train_idx);
  A_test = A(test_idx,:);
  B_test = B(test_idx);
  [X, Y] = inferLowRankImpl(A_train, B_train, tx, rx, lambda, r, mu0, rho, tol_rel, tol_abs, maxiter);
  
  quality = 1 - norm(abs(A_test*X)-B_test)/norm(B_test);

  %----------------------------------------------------------
  % Refine X with additional measurement
  %----------------------------------------------------------
  if quality > 0.6
    X0 = X;
    Y0 = Y;
    [X, Y] = InferADMM(A, B, X0, true, tx, rx, lambda, mu0, rho, tol_rel, tol_abs, maxiter, []);
    similarity = norm(X0'*X)/norm(X0)/norm(X);
    if (similarity < 0.6)
      % rollback to previous X
      X = X0;
      Y = Y0;
    end
  else
    [X, Y] = InferADMM(A, B, X, true, tx, rx, lambda, mu0, rho, tol_rel, tol_abs, maxiter, []);    
  end

  %----------------------------------------------------------
  % scale the solution back
  %----------------------------------------------------------
  X = X * (B_norm / A_norm);
  Y = Y * (B_norm / A_norm);
  
end

function [X, Y, converged] = inferLowRankImpl(A, B, tx, rx, lambda, r, mu0, rho, tol_rel, tol_abs, maxiter)
% 
%   Find X, Y, Z
%
%   minimize:   1/2*|abs(Y) - B|_F^2 + lambda/2*|X|_F^2 + I(Z)   (1)
%
%   subject to:  A*X=Y                                           (2)
%                  X=Z
%
%   where:
%     A: matrix of size T-by-(TX * RX) 
%     X: vector of size (TX * RX)-by-1
%     Y: vector of size T-by-1
%     Z: matrix of size TX-by-RX
%     B: vector of size T-by-1
%
%     abs(Y) is the magnitude vector of complex vector Y
%     B is the measured magnitude of Y
%     I(Z) = 0      if best rank-rz approximation to Z account 
%                   for at least 95% variance of Z
%          = inf    otherwise
%   
%   Input:
%     A:        T-by-(TX * RX) complex matrix
%     B:        T-by-1 RSSI measurements (real numbers)
%     tx:       TX
%     rx:       RX
%     r: upper bound on rank(reshape(Z, tx, rx))
%     tol:      convergence tolerance
%
%   Output:
%     X: [TX*RX]-by-1 complex channel vector
%     Y: T-by-1 complex vector
%     Z: TX-by-RX complex matrix
% 
%   We develop an Alternating Direction Multiplier Method (ADMM)
%   motivated by the following two papers:
%
%   The Augmented Lagrange Multiplier Method for
%   Exact Recovery of Corrupted Low-Rank Matrices
%   http://decision.csl.illinois.edu/~yima/psfile/Lin09-MP.pdf
%
%   Alternating Direction Algorithms for ``$\ell_1$-Problems in 
%   Compressive Sensing
%   http://www.caam.rice.edu/~yzhang/reports/tr0937.pdf
% 
%   Consider the augmented Lagrangian
% 
%      L(X,Y,Z,M,N,mu) =
%        1/2*|abs(Y) - B|_2^2 + lambda/2*|X|_F^2 + I(Z)
%         + <M, A*X-Y>
%         + <N, X-Z>
%         + mu/2 |A*X-Y|_2^2
%         + mu/2 |X-Z|_2^2                             (3)
% 
%    where M is the Lagrangian multipler for constraint A*X=Y
%          N is the Lagrangian multipler for constraint X=Z
% 
%    Alternating minimization of (3) in X, Y, Z and subsequent maximization
%    of M, N via a step of gradient ascent yields the following ADMM:
% 
%         X^{k+1} = argmin_X L(X,Y^{k},Z^{k},M^{k},N^{k},mu)        (4)
%
%         Z^{k+1} = argmin_Z L(X^{k+1},Y^{k},Z,M^{k},N^{k},mu)      (5)
%         Y^{k+1} = argmin_Y L(X^{k+1},Y,Z^{k},M^{k},N^{k},mu)      (6)
%
%         M^{k+1} = M^{k} + mu*(A*X - Y)                            (7)
%         N^{k+1} = N^{k} + mu*(X - Z)                              (8)
% 
%    Next we show how to solve (4), (5) and (6)
% 
% * Find X to minimize L(X,Y,Z,M,N,mu) given Y, Z, M, N, mu:
% 
%   This is simply a least squares problem:
% 
%      minimize mu/2*|AX-Y+M/mu|_2^2 + mu/2*|X-Z+N/mu|_2^2 + lambda/2*|X|_2^2
% 
%   The solution is:
% 
%      x = inv(A'A + I * (1+lambda/mu)) (A'*(Y-M/mu) + Z-N/mu)      (9)
% 
%    where A' is the conjugate transpose of matrix A
% 
% * Find Z to minimize L(X,Y,Z,M,N,mu) given X, Y, M, N, mu:
%
%     minimize I(Z) + mu/2*|X-Z+N/mu|_2^2
%
%   Let E = reshape(X+N/mu, tx, rx).  It becomes:
%   
%     minimize I(Z) + mu/2*|Z-E|_2^2
%
%   Z is simply the best rank-(r) approximation to E, which can be
%   solved via SVD (by keeping only the r largest singular values)
%
% * Find Y to minimize L(X,Y,Z,M,N,mu) given X, Z, M, N, mu:
% 
%     minimize 1/2*|abs(Y)-B|_2^2 + mu/2*|AX-Y+M/mu|_2^2
% 
%   Let C = AX + M/mu. The objective is separable w.r.t. each Y_i:
% 
%     minimize:    1/2 * (abs(Y_i)-B_i)^2 
%               + mu/2 * |Y_i - C_i|^2                             (11)
% 
%   Let R_i = abs(Y_i) be Y_i's magnitude.  It is clear that
%   for any given R_i >= 0, (9) is minimized when Y_i =
%   C_i/||C_i||*R_i. That is, Y_i has the same direction as
%   C_i, but with magnitude r_i.
% 
%   To find R_i, let D_i = ||B_i||.  (9) becomes:
% 
%     minimize:    F(R) =    1/2 * (R-B_i)^2
%                         + mu/2 * (R-D_i)^2                       (12)
%     subject to:  R >= 0
% 
%   F(R) is a degree-2 polynomial of R.  Its minimizer is:
% 
%     R = (B_i + mu*D_i) / (1+mu)
% 
% * Initialization via Over-parameterization.
% 
%   Instead of searching for tx*rx-by-1 vector X, we search for a 
%   tx*rx-by-r matrix instead.  After it converges, we then use
%   the first principle component of X as the initial solution
%   and reset r = 1.
%
  
  [m,n] = size(A);
  [m,t] = size(B);  
  r = min([r m n]);

  if lambda == 0
    U = inv(A'*A + eye(n,n));
    D = [];
  else
    [U,D] = eig(A'*A);
    D = max(0,real(D));
  end

  %----------------------------------------------------------
  % step 1: spectral initialization
  %----------------------------------------------------------
  As = A;
  for i = 1:m
    an = norm(A(i,:));
    if an ~= 0
      As(i,:) = A(i,:) * (B(i) / an);
    end
  end
  AtA = As'*As;
  [V,S] = eig(AtA);
  s2 = max(0,real(diag(S)));
  [s2,idx] = sort(s2,'descend');
  X = bsxfun(@times, V(:,idx(1:r)), sqrt(s2(1:r))');

  %----------------------------------------------------------
  % step 2: over-parameterization
  %----------------------------------------------------------
  scale_by_row = true;
  [X, Y] = InferADMM(A, B, X, scale_by_row, tx, rx, lambda, mu0, rho, tol_rel, tol_abs, maxiter, U, D);
  
  %----------------------------------------------------------
  % step 3: orthonormalize columns of X
  %----------------------------------------------------------
  [Vx,Dx] = eig(X'*X);
  X = X*Vx;
  
  %----------------------------------------------------------
  % step 4: parallel refinement
  %----------------------------------------------------------
  scale_by_row = false;
  [X, Y, converged] = InferADMM(A, B, X, scale_by_row, tx, rx, lambda, mu0, rho, tol_rel, tol_abs, maxiter, U, D);
end

% 
%   Find X, Y
%
%   minimize:   1/2*|abs(Y) - B|_2^2 + I(Z)              (1)
%
%   subject to:  A*X=Y                                   (2)
%                X=Z
%
function [X,Y,converged] = InferADMM(A, B, X0, scale_by_row, tx, rx, lambda, mu0, rho, tol_rel, tol_abs, maxiter, U, D)

  [m,n] = size(A);
  r = size(X0,2);
  
  if isempty(U)
    if lambda == 0
      U = inv(A'*A + eye(n,n));
      D = [];
    else
      [U,D] = eig(A'*A);
      D = max(0,real(D));
    end
  end

  M = zeros(m,r);
  N = zeros(n,r);
  X = X0;
  AX = A*X;
  if scale_by_row
    X = X * (norm(B)/norm(AX,'fro'));
  else
    for j = 1:r
      X(:,j) = X(:,j) * (norm(B)/norm(AX(:,j)));
    end
  end
  AX = A*X;
  Y = normalize_rows(AX, B, scale_by_row);
  Z = ArgMinZ(X, N, 1, tx, rx, m, n);
  AtY = A'*Y;
  
  % Step 2: iteration
  mu = mu0;
  opt_obj = inf;
  converged = false;
  last_res = inf;

  for iter = 1:maxiter

    Y0 = Y;
    Z0 = Z;
    AtY0 = AtY;
    
    % update X
    X = ArgMinX(A, Y, Z, M, N, mu, lambda, U, D);
    AX = A*X;
    
    % update Y
    Y = ArgMinY(AX, B, M, mu, scale_by_row);
    AtY = A'*Y;
    
    % update Z
    Z = ArgMinZ(X, N, mu, tx, rx, m, n);
    
    % update M
    J_M = AX - Y;
    M = M + mu*J_M;

    % update N
    J_N = X - Z;
    N = N + mu*J_N;
    
    % save the best solution so far
    if scale_by_row
      obj = norm(sqrt(sum(abs(AX).^2,2)) - B);
      if obj < opt_obj
        opt_obj = obj;
        opt_X = X;
        opt_Y = Y;
        %disp(['iter: ', int2str(iter), ' r: ', int2str(r), ' obj: ', num2str(opt_obj)]);
      end
    else
      objs = sqrt(sum(bsxfun(@minus, abs(AX), B).^2,1));
      [obj,j] = min(objs);
      if (obj < opt_obj)
        opt_obj = obj;
        opt_X = X(:,j);
        opt_Y = Y(:,j);
        %disp(['iter: ', int2str(iter), ' r: ', int2str(r), ' obj: ', num2str(opt_obj)]);
      end
    end
    
    % convergence test
    res_prim = sqrt(norm(J_M,'fro')^2 + norm(J_N,'fro')^2);
    res_dual = mu*sqrt(norm(AtY - AtY0,'fro')^2 + norm(Z-Z0,'fro')^2);
    res_comb = sqrt(res_prim^2 + norm(Y-Y0,'fro')^2 + norm(Z-Z0,'fro')^2);

    thresh_prim = tol_abs * sqrt((m+n)*r) + tol_rel * sqrt(max(norm(AX,'fro'), norm(Y,'fro'))^2 + max(norm(X,'fro'), norm(Z,'fro'))^2);
    thresh_dual = tol_abs * sqrt(n*r*2) + tol_rel * sqrt(norm(AtY,'fro')^2 + norm(Z,'fro')^2); 
    thresh_comb = tol_abs * sqrt((m+n)*r*2) + tol_rel * sqrt(max(norm(AX,'fro'), norm(Y,'fro'))^2 + max(norm(X,'fro'), norm(Z,'fro'))^2 + norm(Y,'fro')^2 + norm(Z,'fro')^2);

    if (res_prim < thresh_prim && res_dual < thresh_dual) || (res_comb < thresh_comb)
      converged = true;
      break;
    end

    % increase mu whenever combined residual makes too little progress
    % this effectively increase the convexity of the augmented Lagrangian L.
    if res_comb > last_res * 0.9
      mu = mu * rho;
    end
    last_res = res_comb;
  end
  X = opt_X;
  Y = opt_Y;
end


% * Find X to minimize L(X,Y,Z,M,N,mu) given Y, Z, M, N, mu:
% 
%   This is simply a least squares problem:
% 
%   minimize mu/2*|AX-Y+M/mu|_2^2 + mu/2*|X-Z+N/mu|_2^2 + lambda/2*|X|_2^2
% 
%   The solution is:
% 
%      x = inv(A'A + I*(1+lambda/mu)) (A'*(Y-M/mu) + (Z-N/mu))        (9)
% 
%    where A' is the conjugate transpose of matrix A
%
function [X] = ArgMinX(A, Y, Z, M, N, mu, lambda, U, D)
  if lambda == 0
    % U == inv(A'*A+I)
    X = U*(A'*(Y-M/mu) + (Z-N/mu));
  else
    D_inv = 1./(diag(D) + (1 + lambda/mu));
    X = U*bsxfun(@times, D_inv, U'*(A'*(Y-M/mu) + (Z-N/mu)));
  end
end

% * Find Z to minimize L(X,Y,Z,M,N,mu) given X, Y, M, N, mu:
%
%     minimize I(Z) + mu/2*|X-Z+N/mu|_2^2
%
%   Let E = reshape(X+N, tx, rx).  It becomes:
%   
%     minimize I(Z) + mu/2*|Z-E|_2^2
%
%   Let [Ue,Se,Ve] = svd(E).
%
%     Z = Ue * low_rank_approx(Se) * Ve'                          (10)
%
function [Z] = ArgMinZ(X, N, mu, tx, rx, m, n)
  Z = X+N/mu;

  E = reshape(Z, tx, []);

  [U,S] = eig(E*E');
  s2 = max(0,real(diag(S)));
  [s2,idx] = sort(s2,'descend');

  %
  % We enforce a set of low-rank constraints C(r,f): 
  %
  %   first r dimensions capture >= fraction f of variance
  %
  sz = min(rx, tx);
  r0 = ceil(sqrt(sz)*0.5);
  r1 = ceil(sqrt(sz)*0.7);
  r2 = ceil(sqrt(sz));
  r3 = min(sz, ceil(sqrt(sz)*2.0));
  
  f0 = 0.8;
  f1 = 0.9;
  f2 = 0.95;
  f3 = 0.995;

  if m >= n * 3
    r_list = [r3];
    f_list = [f3];
  elseif r1 <= 2
    % sz too small => fallback to simple profile
    r_list = [r2];
    f_list = [f2];
  elseif r0 <= 2
    r_list = [r1 r2 r3];
    f_list = [f1 f2 f3];
  else
    r_list = [r0 r1 r2 r3];
    f_list = [f0 f1 f2 f3];   
  end
  
  %
  % Perform singular value rescaling to satisfy the above low-rank constraints
  %
  s2_scale = ones(size(s2));
  for k = 1:length(r_list)
    r = r_list(k);
    f = f_list(k);
    vr = sum(s2(1:r));
    v = sum(s2);
    if vr < v * f
      scale = min(1, vr / (v - vr) * (1/f-1));
      s2((r+1):end) = s2((r+1):end) * scale;
      s2_scale(idx((r+1):end)) = s2_scale(idx((r+1):end)) * scale;
    end
  end
  
  if (any(s2_scale < 1))
    Z = reshape(bsxfun(@times, U, sqrt(s2_scale)')*U'*E, tx*rx, []);
  end
end

% * Find Y to minimize L(X,Y,M,mu) given X, M, mu:
% 
%     minimize 1/2*|abs(Y)-B|_2^2 + mu/2*|AX-Y+M/mu|_2^2
% 
%   Let C = AX + M/mu. The objective is separable w.r.t. each y_i:
% 
%     minimize:    1/2 * (abs(Y_i)-B_i)^2 
%               + mu/2 * |Y_i - C_i|^2                  (8)
% 
%   Let R_i = abs(Y_i) be Y_i's magnitude.  It is clear that
%   for any given R_i >= 0, (9) is minimized when Y_i =
%   C_i/||C_i||*R_i. That is, Y_i has the same direction as
%   C_i, but with magnitude r_i.
% 
%   To find R_i, let D_i = ||C_i||.  (9) becomes:
% 
%     minimize:    F(R) =    1/2 * (R-B_i)^2
%                         + mu/2 * (R-D_i)^2            (9)
%     subject to:  R >= 0
% 
%   F(R) is a degree-2 polynomial of R.  Its minimizer is:
% 
%     R = (B + mu*D) / (1+mu)
% 
function [Y] = ArgMinY(AX, B, M, mu, scale_by_row)
  Y = AX + M/mu;
  r = size(Y,2);
  if scale_by_row
    D = sqrt(sum(abs(Y).^2,2));
    I = find(D==0);
    if ~isempty(I)
      Y(I,:) = 1/sqrt(r);
      D(I) = 1;
    end
    BD = B./D;
    Y = bsxfun(@times, Y, (BD+mu)/(1+mu));
  else
    D = abs(Y);
    I = find(D==0);
    if ~isempty(I)
      Y(I) = 1;
      D(I) = 1;
    end
    BD = bsxfun(@rdivide, B, D);
    Y = Y.*((BD+mu)/(1+mu));
  end
end

%
% Scale rows of Y.
%
function [Y] = normalize_rows(Y, B, scale_by_row)
  r = size(Y,2);
  if scale_by_row
    D = sqrt(sum(abs(Y).^2,2));
    I = find(D==0);
    if ~isempty(I)
      Y(I,:) = 1/sqrt(r);
      D(I) = 1;
    end
    BD = B./D;
    Y = bsxfun(@times, Y, BD);
  else
    D = abs(Y);
    I = find(D==0);
    if ~isempty(I)
      Y(I) = 1;
      D(I) = 1;
    end
    BD = bsxfun(@rdivide, B, D);
    Y = Y.*BD;
  end
end
