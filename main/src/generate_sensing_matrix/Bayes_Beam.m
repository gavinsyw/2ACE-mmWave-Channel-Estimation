function picked_beams = Bayes_Beam(M,ULA,cb_train,num_sel)

Nt = ULA.Nt;
Nr = ULA.Nr;

% randomlh pick some beams first, to narrow search space
% candidate_size = 100*M;
candidate_size = 40000;
candidate_set = randi(num_sel, candidate_size, 1);
cb_candidate = cb_train(candidate_set, :);

[picked_candidates, ~] = bayesAopt_complex(cb_candidate, M, 'K', eye(Nt*Nr));

picked_beams = candidate_set(picked_candidates);

end