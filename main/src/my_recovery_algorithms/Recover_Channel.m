function H = Recover_Channel(measurements, beams, Method, ULA, s, AD)

if Method.ADMM
    disp('Estimating with ADMM');
    version = 0;
    [H.ADMM, ~, ~] = ADMM_v2(measurements, beams, ULA.Nt, ULA.Nr, version);
end

if Method.ADMMLowRankV1
    disp('Estimating with ADMMLowRankV1');
    version = 1;
    [H.ADMMLowRankV1, ~, ~] = ADMM_v2(measurements, beams, ULA.Nt, ULA.Nr, version);
end

if Method.ADMMLowRankV2
    disp('Estimating with ADMMLowRankV2');
    version = 2;
    [H.ADMMLowRankV2, ~, ~] = ADMM_v2(measurements, beams, ULA.Nt, ULA.Nr, version);
end

if Method.ADMMLowRankV3
    disp('Estimating with ADMMLowRankV3');
    version = 3;
    [H.ADMMLowRankV3, ~, ~] = ADMM_v2(measurements, beams, ULA.Nt, ULA.Nr, version);
end

if Method.ADMMLowRankV4
    disp('Estimating with ADMMLowRankV4');
    version = 4;
    [H.ADMMLowRankV4, ~, ~] = ADMM_v2(measurements, beams, ULA.Nt, ULA.Nr, version);
end

if Method.PhaseLift
    disp('Estimating with PhaseLift');
    H.PhaseLift = MyPhaseLift((measurements./2e5).^2.*1e10, beams)./sqrt(1e10).*2e5;
end

if Method.PLOMP || Method.PLGAMP
    disp('Estimating with PLOMP and PLGAMP');
    [plomp_signal, plgamp_signal] = ...
        My_TwoStage_Recovery((measurements./2e5).^2.*1e10, beams*AD, s, 1, 0, Method);
    % recover channel from signal
    H.PLOMP = AD * plomp_signal./sqrt(1e10).*2e5;
    H.PLGAMP = AD * plgamp_signal./sqrt(1e10).*2e5;
end

end
