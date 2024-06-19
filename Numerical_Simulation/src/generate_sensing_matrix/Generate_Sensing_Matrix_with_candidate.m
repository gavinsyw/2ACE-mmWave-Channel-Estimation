function picked_beams = Generate_Sensing_Matrix_with_candidate(Beampattern_Mode,M,ULA,cb_train,rss_train)
%% Parameter fetching
% number of selections and size of FW
[num_sel, ~] = size(cb_train);
% AD = Sparse_Channel_Representation.AD;
% [~, AD_len] = size(AD);

%% Generate beam patterns 
switch Beampattern_Mode
    case 'Random_Phase_State' % randomly pick some beams in the given codebook
%         picked_beams = randi(num_sel, M, 1);
        picked_beams = 1:M;
% 
%     case 'Directional_Beam'                     % directional beam pattern with uniform gain in spatial domain
%         Rank_Eliminated = 0;
%         [F, W_single] = Directional_Beam(Mt,Mr,ULA,H,Rank_Eliminated);
%         W = zeros(U, Nr, Mr);
%         FW = zeros(U, Mt*Mr, Nt*Nr);
%         A = zeros(U, Mt*Mr, AD_len);
%         for i = 1:U
%             W(i,:,:) = W_single;
%             FW(i,:,:) = kron(transpose(F),W_single');
%             A(i,:,:) = squeeze(FW(i,:,:))*AD;
%         end
% 
%     case 'Directional_Beam_Angular'             % directional beam pattern with uniform gain in angular domain
%         [F, W_single] = Directional_Beam_Angular(Mt,Mr,ULA,H);   
%         W = zeros(U, Nr, Mr);
%         FW = zeros(U, Mt*Mr, Nt*Nr);
%         A = zeros(U, Mt*Mr, AD_len);
%         for i = 1:U
%             W(i,:,:) = W_single;
%             FW(i,:,:) = kron(transpose(F),W_single');
%             A(i,:,:) = squeeze(FW(i,:,:))*AD;
%         end

    case 'Bayes_Beam'
        picked_beams = Bayes_Beam(M,ULA,cb_train,num_sel);

    otherwise
        ;
end


end