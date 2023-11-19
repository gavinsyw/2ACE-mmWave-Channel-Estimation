function Evaluation_result = Evaluate_H(recoveredSig,L,H,Sparse_Channel_Representation,ULA,SNR,Mtr,noise_power,U)
    % load parameters
    Nt = ULA.Nt;
    Nr = ULA.Nr;
    
    % estimate normalize error
    vecH_Estimated = recoveredSig;
    phaseFac = exp( 1i* angle( (vecH_Estimated'*transpose(H.vecH(1,:)))/(H.vecH(1,:)*H.vecH(1,:)') ) );
    vecH_Estimated = vecH_Estimated*phaseFac;
    H_estimated = reshape(vecH_Estimated,Nr,Nt);
    H_real = reshape(H.H_Matrix(1,:,:),Nr,Nt);
    MSE_H = norm(H_estimated-H_real,'fro')^2/norm(H_real,'fro')^2;
    
    X_gt = transpose(H.vecH);
    X = vecH_Estimated;
    MSE_H = norm(X_gt - (X'*X_gt)/(X'*X)*X,'fro')^2/norm(X_gt,'fro')^2;
    
    disp(MSE_H);
    
    Evaluation_result=[0;...
                   0;...
                   0;...
                   MSE_H];
end