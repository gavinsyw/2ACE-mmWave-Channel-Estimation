function H = Recover_Channel(measurements, beams, Method, ULA, s, AD)

if Method.ADMM
    version = 0;
    [H.ADMM, ~, ~] = ADMM_v2(measurements, beams, ULA.Nt, ULA.Nr, version);
end

if Method.ADMMLowRankV1
    version = 1;
    [H.ADMMLowRankV1, ~, ~] = ADMM_v2(measurements, beams, ULA.Nt, ULA.Nr, version);
end

if Method.ADMMLowRankV2
    version = 2;
    [H.ADMMLowRankV2, ~, ~] = ADMM_v2(measurements, beams, ULA.Nt, ULA.Nr, version);
end

if Method.ADMMLowRankV4
    version = 3;
    [H.ADMMLowRankV4, ~, ~] = ADMM_v2(measurements, beams, ULA.Nt, ULA.Nr, version);
end

if Method.PLOMP || Method.PLGAMP
    [plomp_signal, plgamp_signal] = ...
        My_TwoStage_Recovery(measurements.^2, beams*AD, s, 1, 0, Method);
    % recover channel from signal
    H.PLOMP = AD * plomp_signal;
    H.PLGAMP = AD * plgamp_signal;
end

if Method.PhaseLift
    pl_signal = MyPhaseLift(measurements.^2, beams);
    % recover channel from signal
    H.PhaseLift = pl_signal;
end

end