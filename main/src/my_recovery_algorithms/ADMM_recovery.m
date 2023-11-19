function [recovered_channel, recovered_signal, multiplier] = ADMM_recovery(measurements, F, lam, learning_rate)

c = measurements.^2;
A = transpose(F);
mu = learning_rate;

[n_tries, n_tx] = size(A);

% initialization
x = zeros(n_tx, 1);     % channel vector
y = zeros(n_tries, 1);  % recovered signal
m = zeros(n_tries, 1);  % Lagrangian multiplier

for i = 1:n_tries
    % fixed amplitude, random phase
    y(i) = sqrt(c(i)) * exp(1j*rand()*2*pi);
end

% Iterate for fixed number of times. can be further improved to detect
% convergence.
for k = 1:100
    % solve x
    x = (A'*A+lam*eye(n_tx))\A'*y;
    
    % solve y
    b = A*x;
    ang_y = angle(b);
    e = abs(b);
    d = c - m/mu;
    amp_y = zeros(n_tries, 1);
    for i = 1:n_tries
        % find amplitude of y using equation solver
        syms r;
        s = solve(2*mu*r^3+(1-2*mu*d(i))*r-e(i)==0, r, 'Real', true);
        amp_y(i) = abs(s(1));
    end
    y = amp_y .* exp(1j*ang_y);
    
    % dual update
    m = m + mu*(abs(y).^2-c);
end

recovered_channel = x;
recovered_signal = y;
multiplier = m;

end