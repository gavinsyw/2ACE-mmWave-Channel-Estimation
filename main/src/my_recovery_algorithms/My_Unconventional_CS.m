function recovered_H = My_Unconventional_CS(measurements, F)
    b = measurements;
    A = transpose(F);
    [~, w] = size(A);
    
    % singular value decomposition
    [U, S, ~] = eig(A'*A);
    c = U'*A'*b;
    s = diag(S);
    
    fun_lam = @(lam)(norm(c./(s+lam))-1)^2;
    
    lambda = fmincon(fun_lam, 0.5, [], [], [], [], 0, 1);
    
    recovered_H = (A'*A+lambda*eye(w))^(-1)*A'*b;
end