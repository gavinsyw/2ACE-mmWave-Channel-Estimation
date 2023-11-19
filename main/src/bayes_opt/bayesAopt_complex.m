function [rowlist,Acriteria] = bayesAopt_complex(fxcand,nrows,varargin)
%BAYESAOPT Bayesian A-optimal design from candidates using row exchanges.
%
%   The high level structure of the code is similar to bayesDopt
%   thought the detail is different
%
%   RLIST = BAYESAOPT(CAND,NROWS) uses a row-exchange algorithm to 
%   select a Bayesian A-optimal design from the candidate set CAND.  
%   CAND is an N-by-P matrix containing the values of P model terms at each of 
%   N points.  NROWS is the desired number of rows in the design.  
%   RLIST is a vector of length NROWS listing the selected rows.
%
%   The BAYESAOPT function selects a starting design X at random, and
%   uses a row-exchange algorithm to iteratively replace rows of X by
%   rows of CAND in an attempt to reduce the Bayesian A-criteria
%   trace(A*inv(X'*X+K)).
%
%   RLIST = BAYESAOPT(CAND,NROWS,'PARAM1',VALUE1,'PARAM2',VALUE2,...)
%   provides more control over the design generation through a set of
%   parameter/value pairs.  Valid parameters are the following:
%
%      Parameter    Value
%      'init'       Initial design as an NROWS-by-P matrix (default
%                   is a random subset of the rows of CAND).
%      'maxiter'    Maximum number of iterations (default = 10).
%      'c'          C in Bayesian A-criterion: trace(C'*C*inv(X'*X+K))
%                   (default: I)
%      'k'          K in Bayesian A-criterion: trace(C'*C*inv(X'*X+K))
%                   (default: I/10)
%      'a'          C'*C (if 'a' is not empty, will disregard 'c')
%      'fixed'      fixed set of rows that cannot change (size(fixed,2)=size(X,2))
%
%
%   Example:  generate a Bayesian D-optimal design when there is a restriction
%   on the candidate set, so the ROWEXCH function isn't appropriate.
%      F = (fullfact([5 5 5])-1)/4;   % factor settings in unit cube
%      T = sum(F,2)<=1.51;            % find rows matching a restriction
%      F = F(T,:);                    % take only those rows
%      C = [ones(size(F,1),1) F F.^2];% compute model terms including
%                                     % a constant and all squared terms
%      R = bayesAopt(C,12);           % find a D-optimal 12-point subset
%      X = F(R,:);                    % get factor settings
%
%   See also CANDEXCH, CANDGEN, ROWEXCH, CORDEXCH, X2FX.
%
%   Technical Details:  
%
%   A-optimality is useful for dealing with quadratic loss function in a linear
%   normal system y = X*theta.  
%
%   Specifically, in the Bayesian A-criteria, A is a symmetric non-negative 
%   definite matrix that is used to weight different estimation problem where 
%   the interest may be in estimating individual components of theta or linear 
%   combination of theta.  The corresponding utility function is simply
%
%         E[ (theta-hat(theta))' A (theta-hat(theta)) ] = cov( C*theta )
%
%   where C'C = A.
%
%   where hat(theta) is the estimated value of theta. For example, if we want
%   to minimize the RMSE of X*theta (i.e. e2e performance), we should set A to
%   (X'*X).
%
%   K  corresponds to R in the standard Bayesian A-criteria, where
%   theta is assumed to have covariance matrix sigma^2 K^{-1} that is, 
%   sigma^{-2} K is the prior precision matrix.
%
% Additional Notes: (XXX: VERY IMPORTANT)
%
%
% (1) an alternative way to understand the R (i.e. K)
%
%   Intuitively, given y = X*theta, let Sigma = cov(theta) and SS'=Sigma,
%   we can then transform the problem to y = X*S*(inv(S)*theta), where
%   inv(S)*theta has identity covariance.  The metric then becomes
%
%          trace( A*Sigma*inv((X*S)'(X*S) + I/n) )
%
%   which is equivalent to 
%
%          trace( A*inv(X'*X + Sigma^-1/n) ) = trace( A*inv(X'X + R/n) )
%
% (2) how to deal with difference variance of y-X*theta?
%
%   Finally, the above analysis assumes that for y = X*theta + epsilon,
%   epsilon has identity covariance.  If this is not the case, then we
%   need to scale the rows of X & A properly s.t. epsilon has identity
%   covariance structure.
%
% (3) Prior distribution (theta0, Sigma)
%
%   For the actual analysis, we first translate the problem into 
%   y = X*(theta-theta0) + theta0, and try to solve y-theta0.
%
% References:
%
%   Bayesian Experimental Design: A Review (1995) 
%   Kathryn Chaloner, Isabella Verdinelli
%   Statist. Sci.
%   http://citeseer.ist.psu.edu/chaloner95bayesian.html

%   Copyright 1993-2002 The MathWorks, Inc. 
%   $Revision: 1.2 $  $Date: 2021/03/08 04:30:16 $

rowlist = zeros(nrows,1);
maxiter = 2;
iter = 0;
madeswitch = 1;
acutoff = -sqrt(eps);
A = [];
fixed = fxcand([],:);

okargs = {'maxiter' 'init' 'c' 'k' 'a' 'fixed'};
defaults = {5 [] speye(size(fxcand,2)) speye(size(fxcand,2))/1000 [] fxcand([],:)};
[emsg,maxiter,xinit,C,K,A,fixed] = statgetargs(okargs,defaults,varargin{:});

C = sparse(C);
if length(A) == 0
  A = C'*C;
end
error(emsg);

if ~isnumeric(maxiter) | prod(size(maxiter))~=1 | maxiter<1
   error('Value of the ''maxiter'' parameter must be a number 1 or greater.');
end

if isempty(xinit)
   X = fxcand(unidrnd(size(fxcand,1),nrows,1),:);
else
   X = xinit;
end

% Some simple matrix algebra 
% (see http://mcs.une.edu.au/~nkn/pdf/review.pdf for details)  
%
% NOTE: http://mcs.une.edu.au/~nkn/pdf/review.pdf has typos on the sign of 
% w and w_i.  What I have below is correct.
%
% Given M = X'*X + K
% 
% If x' is a row vector to be removed from X, then we have
%
%    | M-xx' |   = |M| (1 - x'*inv(M)*x)              (1) // |M| = det(A)
%
%    inv(M-xx')  = inv(M) + wuu'                      (2)
%    where u = inv(M)*x and w = 1/(1-x'inv(M)x) = 1/(1-x'u)
%
% Let N = M - xx', if x_i' is a row vector to be added to the new X,
% we have
%
%    | N + x_i x_i' | = |N| (1 + x_i' inv(N) x_i)     (3)
%
%    inv(N + x_i x_i') = inv(N) - w_i u_i u_i'        (4)
%    where u_i = inv(N)*x_i and w_i = 1/(1 + x_i' inv(N) x_i) 
%    = 1/(1 + x_i' u_i)
%
% Putting them together, we have
%
%   | M - xx' + x_i x_i' | = |M| (1 - x'inv(M)x) (1 + x_i' inv(N) x_i) (5)
%
%                          = |M| {1 + Delta(x_i, x)}                   (6)
%
%   we refer to Delta(x_i, x) as Fedorov's delta function.  This is useful 
%   for computing D-optimal designs.
%
%  For A-optimality, we have
%
%    trace(A*inv(M - xx')) = trace(A*inv(M)) + trace(wAuu')
%                          = trace(A*inv(M)) + w sum((A*u).*conj(u))
%
%    trace(A*inv(N + x_i x_i')) = trace(A*inv(N)) - w_i sum((A*u_i).*conj(u_i))
%

fxcand    = sparse(fxcand);
K         = sparse(K);
A         = sparse(A);
X         = sparse(X);
fixed     = sparse(fixed);

Minv      = inv(X'*X+K+fixed'*fixed);
MinvF     = Minv*fxcand';
AUUtrace  = sum((A*MinvF).*conj(MinvF));
Acriteria = abs(trace(A*Minv));

FMinvF    = abs(sum(conj(fxcand)'.*MinvF));

while madeswitch > 0 & iter < maxiter
  madeswitch = 0;
  iter = iter + 1;
  
  for row = 1:nrows

%     [iter, row]
%     if (mod(row,10) == 0) 
%       fprintf('bayesAopt; iter %d row %d\n', iter, row);
%     end

    x       = X(row,:)';
    u       = Minv*x;
    uf      = (fxcand*u)';

    w       = real(1/(1-x'*u));

    wu      = (w*u);
    wuf     = (w*uf);
    Au      = (A*u);

    %Ninv    = Minv   + u*wu';       % Ninv   = inv(M-xx')     % delayed
    %NinvF   = MinvF  + u*wuf;       % NinvF  = Ninv*fxcand'    
    %ANinvF  = AMinvF + Au*wuf;      % ANinvF = A*Ninv*fxcand'
    FNinvF  = FMinvF + w*abs(uf).^2;

    %U       = NinvF;
    %W       = 1./(1+sum(fxcand'.*U));
    W       = 1./(1+FNinvF);

    AUUtrNew = AUUtrace + 2*real((Au'*Minv*fxcand').*conj(wuf)) + (u'*Au)*abs(wuf).^2;

    %ad     =  w*sum((AMinv*x).*u) - W.*sum(ANinvF.*U);
    %ad      =  sum(Au.*wu) - W.*sum(ANinvF.*U);
    ad      = real(sum(Au.*conj(wu)) - W.*AUUtrNew);
    [a,idx] = min(ad);

    % Switch rows if the maximum change is greater than 1.
    if (a < acutoff) | (rowlist(row) == 0)
      madeswitch = 1;
      Acriteria = Acriteria + a;
      rowlist(row) = idx;
      X(row,:) = fxcand(idx,:);      
      
      % first apply delayed update of Ninv and ANinv
      Ninv    = Minv  + u*wu';       % Ninv   = inv(M-xx')
      NinvF  = Ninv*fxcand';
      
      % then update Minv, AMinv, MinvF, AMinvF      
      x       = X(row,:)';
      u       = Ninv*x;
      w       = W(idx);

      uf      = (fxcand*u)';
      
      wu      = (w*u);
      wuf     = (w*uf);
      Au      = (A*u);
      
      Minv    = Ninv - u*wu';       % Minv   = inv(N+xx')
      %MinvF   = NinvF  - u*wuf;       % MinvF  = Minv*fxcand'
      %AMinvF  = ANinvF - Au*wuf;     % AMinvF = A*Minv*fxcand'
      FMinvF  = FNinvF - w*abs(uf).^2;
      
      AUUtrace = AUUtrNew - 2*real((Au'*Ninv*fxcand').*conj(wuf)) + (u'*Au)*abs(wuf).^2;
    end
  end
end
