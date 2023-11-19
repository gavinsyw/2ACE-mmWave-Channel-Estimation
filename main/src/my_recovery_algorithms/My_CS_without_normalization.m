function recovered_H = My_CS_without_normalization(measurements, F)
    b = measurements;
    A = transpose(F);
    
    recovered_H = (A'*A)^(-1)*A'*b;
end