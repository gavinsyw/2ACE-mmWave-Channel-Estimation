function [X, Y, quality] = inferMinL2(A, B, lambda, r, tol_rel, tol_abs, maxiter)
  if nargin < 3, lambda = 0; end
  if nargin < 4, r = 20; end
  if nargin < 5, tol_rel = 1e-4; end
  if nargin < 6, tol_abs = 1e-8; end
  if nargin < 7, maxiter = 500; end
  
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
  train_idx = randsample(m, ceil(m*0.95));
  test_idx = setdiff(1:m, train_idx);
  A_train = A(train_idx,:);
  B_train = B(train_idx);
  A_test = A(test_idx,:);
  B_test = B(test_idx);
  [X, Y] = inferMinL2Impl(A_train, B_train, lambda, r, tol_rel, tol_abs, maxiter);
  
  quality = 1 - norm(abs(A_test*X)-B_test)/norm(B_test);

  %----------------------------------------------------------
  % Refine X with additional measurement
  %----------------------------------------------------------
  if quality > 0.6
    X0 = X;
    Y0 = Y;
    [X, Y] = InferADMM(A, B, X0, true, lambda, tol_rel, tol_abs, maxiter, []);
    similarity = norm(X0'*X)/norm(X0)/norm(X);
    if (similarity < 0.6)
      % rollback to previous X
      X = X0;
      Y = Y0;
    end
  end

  %----------------------------------------------------------
  % scale the solution back
  %----------------------------------------------------------
  X = X * (B_norm / A_norm);
  Y = Y * (B_norm / A_norm);

end

function [X, Y] = inferMinL2Impl(A, B, lambda, r, tol_rel, tol_abs, maxiter)
% 
%   Find X and Y
%
%   minimize:    1/2|mag(Y) - B|_2^2 + lambda/2*|X|_F^2  (1)
%   subject to:  A*X=Y                                   (2)
% 
%   where:
%     A: matrix of size T-by-(TX * RX) 
%     X: vector of size (TX * RX)-by-1
%     Y: vector of size T-by-1
%     B: vector of size T-by-1
%
%     mag(Y(i,:)) = norm(Y(i,:)) is the magnitude of Y(i,:)
%     B is the measured magnitude of Y;
%   
%   Input:
%     A: T-by-(TX * RX) complex matrix
%     B: T-by-1 RSSI measurements (real numbers)
%
%   Output:
%     X: [TX*RX]-by-1 complex channel vector
%     Y: T-by-1 complex vector
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
%      L(X,Y,M,mu) =
%           1/2*|mag(Y) - B|_2^2 + lambda/2*|X|_F^2
%         + <M, A*X-Y> 
%         + mu/2 |A*X-Y|_2^2                              (3)
% 
%    where M is the Lagrangian multipler for constraint A*X=Y
% 
%    Alternating minimization of (3) in X, Y and subsequent maximization
%    of M via a step of gradient ascent yields the following ADMM:
% 
%         X^{k+1} = argmin_X L(X,Y^{k},M^{k},mu)         (4)
%         Y^{k+1} = argmin_Y L(X^{k+1},Y,M^{k},mu)       (5)
%         M^{k+1} = M^{k} + mu*(A*X - Y)                 (6)
% 
%    Next we show how to solve (4) and (5).
% 
% * Find X to minimize L(X,Y,M,mu) given Y, M, mu:
% 
%   This is simply a least squares problem:
% 
%      minimize mu/2*|AX-Y+M/mu|_2^2 + lambda/2*|X|_F^2
% 
%   The solution is:
% 
%      x = inv(A'A + lambda/mu) (A'*(Y-M/mu))            (7)
% 
%    where A' is the conjugate transpose of matrix A
% 
% * Find Y to minimize L(X,Y,M,mu) given X, M, mu:
% 
%     minimize 1/2*|mag(Y)-B|_2^2 + mu/2*|AX-Y+M/mu|_2^2
% 
%   Let C = AX + M/mu. The objective is separable w.r.t. each y_i:
% 
%     minimize:    1/2 * (mag(Y_i)-B_i)^2 
%               + mu/2 * |Y_i - C_i|^2                  (8)
% 
%   Let R_i = mag(Y_i) be Y_i's magnitude.  It is clear that
%   for any given R_i >= 0, (9) is minimized when Y_i =
%   C_i/||C_i||*R_i. That is, Y_i has the same direction as
%   C_i, but with magnitude r_i.
% 
%   To find R_i, let D_i = ||B_i||.  (9) becomes:
% 
%     minimize:    F(R) =    1/2 * (R-B_i)^2
%                         + mu/2 * (R-D_i)^2            (9)
%     subject to:  R >= 0
% 
%   F(R) is a degree-2 polynomial of R.  Its minimizer is:
% 
%     R = (B_i + mu*D_i) / (1+mu)
% 
% * Spectral Initialization
% 
%   We find the largest eigenvector of sum_i B(i)^2 A(i,:)' * A(i,:)
%   using power method.
%

  [m,n] = size(A);
  [m,t] = size(B);

  if lambda == 0
    U = pinv(A);
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
  if sum(s2(1:r)) >= sum(s2) * 0.9
    r = max(min(find(cumsum(s2) >= sum(s2) * 0.9)), 3);
    r = min([r m n]);
  end
  X = bsxfun(@times, V(:,idx(1:r)), sqrt(s2(1:r))');
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
  if sum(s2(1:r)) >= sum(s2) * 0.9
    r = max(min(find(cumsum(s2) >= sum(s2) * 0.9)), 3);
    r = min([r m n]);
  end
  X = bsxfun(@times, V(:,idx(1:r)), sqrt(s2(1:r))');
  
  %----------------------------------------------------------
  % step 2: over-parameterization
  %----------------------------------------------------------
  scale_by_row = true;
  [X, Y] = InferADMM(A, B, X, scale_by_row, lambda, tol_rel, tol_abs, maxiter, U, D);
  
  %----------------------------------------------------------
  % step 3: orthonormalize columns of X
  %----------------------------------------------------------
  [Vx,Dx] = eig(X'*X);
  X = X*Vx;

  %----------------------------------------------------------
  % step 4: parallel refinement
  %----------------------------------------------------------
  scale_by_row = false;
  [X, Y, converged] = InferADMM(A, B, X, scale_by_row, lambda, tol_rel, tol_abs, maxiter, U, D);

end

function [X, Y, converged] = InferADMM(A, B, X0, scale_by_row, lambda, tol_rel, tol_abs, maxiter, U, D)

  [m,n] = size(A);
  r = size(X0,2);
  
  if isempty(U)
    if lambda == 0
      U = pinv(A);
      D = [];
    else
      [U,D] = eig(A'*A);
      D = max(0,real(D));
    end
  end
  
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
  AtY = A'*Y;

  M = zeros(m,r);

  % Step 2: iteration
  mu = 0.001;
  rho = 1.03;
  opt_obj = inf;
  last_res = inf;

  converged = false;
  for iter = 1:maxiter

    Y0 = Y;
    AtY0 = AtY;
    
    % update X
    X = ArgMinX(A, Y, M, mu, lambda, U, D);
    AX = A*X;
    
    % update Y
    Y = ArgMinY(AX, B, M, mu, scale_by_row);
    AtY = A'*Y;
    
    % update M
    J_M = AX - Y;
    M = M + mu*J_M;

    % save the best solution so far
    if scale_by_row
      obj = norm(sqrt(sum(abs(AX).^2,2)) - B);
      if obj < opt_obj
        opt_obj = obj;
        opt_X = X;
        opt_Y = Y;
        % disp(['iter: ', int2str(iter), ' r: ', int2str(r), ' obj: ', num2str(opt_obj)]);
      end
    else
      objs = sqrt(sum((abs(AX) - repmat(B,1,r)).^2,1));
      [obj,j] = min(objs);
      if (obj < opt_obj)
        opt_obj = obj;
        opt_X = X(:,j);
        opt_Y = Y(:,j);
        % disp(['iter: ', int2str(iter), ' r: ', int2str(r), ' obj: ', num2str(opt_obj)]);
      end
    end

    % convergence test
    res_prim = norm(J_M,'fro');
    res_dual = mu*norm(AtY - AtY0,'fro');
    res_comb = sqrt(res_prim^2 + norm(Y-Y0,'fro')^2);

    thresh_prim = tol_abs * sqrt(m*r) + tol_rel * max(norm(AX,'fro'), norm(Y,'fro'));
    thresh_dual = tol_abs * sqrt(n*r) + tol_rel * norm(AtY,'fro');
    thresh_comb = tol_abs * sqrt(m*r*2) + tol_rel * sqrt(max(norm(AX,'fro'), norm(Y,'fro'))^2 + norm(Y,'fro')^2);

    if (res_prim < thresh_prim && res_dual < thresh_dual) || (res_comb < thresh_comb)
      converged = true;
      break
    end    
    
    % increase mu whenever combined residual makes too little progress
    if res_comb > last_res * 0.9
      mu = mu * rho;
    end
    last_res = res_comb;
  end
  
  X = opt_X;
  Y = opt_Y;
end

% * Find X to minimize L(X,Y,M,mu) given Y, M, mu:
%   This is simply a least squares problem:
% 
%      minimize mu/2*|AX-Y+M/mu|_2^2 + lambda/2*|X|_2^2
% 
%   The solution is:
% 
%      x = inv(A'A+I*lambda/mu)*(A'*(Y-M/mu))
%
function [X] = ArgMinX(A, Y, M, mu, lambda, U, D)
  if lambda == 0
    % U == pinv(A)
    X = U*(Y-M/mu);
  else
    D_inv = 1./(diag(D) + lambda/mu);
    X = U*bsxfun(@times, D_inv, U'*(A'*(Y-M/mu)));
  end
end

% * Find Y to minimize L(X,Y,M,mu) given X, M, mu:
% 
%     minimize 1/2*|mag(Y)-B|_2^2 + mu/2*|AX-Y+M/mu|_2^2
% 
%   Let C = AX + M/mu. The objective is separable w.r.t. each y_i:
% 
%     minimize:    1/2 * (mag(Y_i)-B_i)^2 
%               + mu/2 * |Y_i - C_i|^2                  (8)
% 
%   Let R_i = mag(Y_i) be Y_i's magnitude.  It is clear that
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
