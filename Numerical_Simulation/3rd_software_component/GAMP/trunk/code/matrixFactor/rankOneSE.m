% rankOneSE:  State-evolution analysis of rank-one factorization
%
% Performs the SE analysis assuming MMSE estimates matched to the true
% prior.  

% Set path 
addpath('../main');
addpath('../stateEvo');

% Indices for the distributions on U and V
U_DIST = 1;
V_DIST = 2;
ndist = 2;

% Distributions for U and V are set in the distType vector
GAUSS_DIST = 1; % Gaussian
CONST_DIST = 2;
EXP_DIST = 3;
distType = [GAUSS_DIST EXP_DIST];
sparseRat = [1 0.1]';   % Sparsity ratio 

% Other parameters
m = 1000;       % Dimension of u
n = 500;        % Dimemnsion of v
beta = n/m;     % Measurement ratio
snrTest = linspace(-5,15,21)';      % Set of SNRs to compute
nsnr = length(snrTest);             % num SNR points
savedat = 1;    % Flag indicating if results are to be saved at the end

% Generate the EstimInAvg classes for the U and V distributions.  This is
% used for the SE analysis
estim0 = cell(ndist,1);
estim = cell(ndist,1);
estimAvg = cell(ndist,1);
for idist = 1:ndist
    
    % Gaussian distribution
    if (distType(idist) == GAUSS_DIST)
        if (idist == 1) % u
            mean0 = 0;
            var0 = 1;
        else
            mean0 = 0.1; % v
            var0 = 1;
        end
        estim0{idist} = AwgnEstimIn(mean0, var0);
    elseif (distType(idist) == CONST_DIST)
        x = [0 1]';
        px = [0.9 0.1]';
        estim0{idist} = DisScaEstim(x,px);
        
    elseif (distType(idist) == EXP_DIST)
        % Exponential distribution
        nx = 100;
        x = linspace(1/nx,2,nx)';
        px = exp(-x);
        px = px/sum(px);
        estim0{idist} = DisScaEstim(x,px);
    end
    
    % Sparsify the distribution if necessary
    if (sparseRat(idist)<1)
        x = [0; x];
        px = [ 1-sparseRat(idist); sparseRat(idist)*px];
        estim{idist} = SparseScaEstim( estim0{idist}, sparseRat(idist) );
    else
        estim{idist} = estim0{idist};
    end
    
    % Generate the EstimInAvg class which computes the 
    % average MSE for the SE recursions.
    if (distType(idist) == GAUSS_DIST)
        % For the Gaussian estimator, use the pre-built class
        estimAvg{idist} = AwgnEstimInAvg(var0);
    else
        % For the other distributions, compute via Monte Carlo sampling
        nw = 100;
        estimAvg{idist} = MCEstimInAvg( estim{idist}, x, px, nw);
    end
    
    
end

% Set distributions
estimu = estim{1};
estimv = estim{2};

% Get initial variances
[umean0,uvar0] = estimu.estimInit();
[vmean0,vvar0] = estimv.estimInit();
vsq0 = vmean0^2+vvar0;
usq0 = umean0^2+uvar0;

% Initialize vectors
nit = 10;
corru = zeros(nit,nsnr);
corrv = zeros(nit+1,nsnr);
for isnr = 1:nsnr
    
    % Set noise level
    snr = snrTest(isnr);
    wvar = usq0*vsq0*10^(-0.1*snr);
    
    % Get average estimators
    corrv(1,isnr) = vmean0^2/vsq0;
    for it = 1:nit
        
        snru = beta*vsq0/wvar*corrv(it,isnr);
        corru(it,isnr) = 1-estimAvg{1}.avgMSE(1/snru)/usq0;
        
        snrv = usq0/wvar*corru(it,isnr);
        corrv(it+1,isnr) = 1-estimAvg{2}.avgMSE(1/snrv)/vsq0;
        
    end
    
    fprintf(1,'%d snr=%f corrv=%f\n', isnr, snr, corrv(end,isnr) );
end

plot(snrTest, corrv(end,:)');
grid on;

if savedat
    save data/rankOneSE snrTest corrv corru;
end

