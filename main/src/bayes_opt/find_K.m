function K = find_K(H, H_mtx, n, ULA)
    lambda = ULA.lambda;
    d = ULA.d;
    Nt = ULA.Nt;
    Nr = ULA.Nr;
    vecK = zeros(n, 1);
    AoD_range = H.AoD_range(1):(H.AoD_range(2)-H.AoD_range(1))/n:H.AoD_range(2);
    for i = 1:n
        w = transpose(exp(-1i*2*pi*d/lambda*sind(AoD_range(i)).*(0:1:Nt-1)))/sqrt(Nt);
        vecK(i) = sqrt(abs(H_mtx*w)).^-1;
    end
    K = diag(vecK);
end