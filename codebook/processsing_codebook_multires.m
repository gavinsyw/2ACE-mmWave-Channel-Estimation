%###########################################################################
%Project 2ACE
%PI: Prof. Lili Qiu @ UT Austin/MSRA
%Participants: Yiwen Song @ CMU, Changhan Ge @ UT Austin, Yin Zhang
%Copyright @ The University of Texas at Austin, 2023
%Partially inherited from Teng Wei/Song Wang/Jingqi Huang/Xinyu Zhang @ UCSD
%###########################################################################

clear
clc

%% Rx codebook
fid = fopen('./codebook_brd/multires_16ant_rx/multires_16_rx.txt');
tline = fgetl(fid);
rx_codebook_phase = [];
while ischar(tline)
    if tline == -1
        break;
    end
    rx_codebook_phase = [rx_codebook_phase;sscanf(tline,'%1d' )'];
    tline = fgetl(fid);
end

%% Tx codebook
fid = fopen('./codebook_brd/multires_16ant_tx/multires_16_tx.txt');
tline = fgetl(fid);
tx_codebook_phase = [];
counter = 0;
while ischar(tline)
    counter = counter+1;
    if tline == -1
        break;
    end
    cb = sscanf(tline,'%1d' )';
    if ~any(cb) && (mod(counter,64)==1 || mod(counter,64)==2)
        tline = fgetl(fid);
        continue;
    end
    tx_codebook_phase = [tx_codebook_phase;cb];
    tline = fgetl(fid);
end

codebook_amp = [1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0];
id = [1 2 3 4 5 6 7 8 17 18 19 20 21 22 23 24];

tx_codebook = codebook_amp.* exp(1j.*(pi./2.*tx_codebook_phase));
tx_codebook = tx_codebook(:,id);
rx_codebook = codebook_amp.* exp(1j.*(pi./2.*rx_codebook_phase));
rx_codebook = rx_codebook(:,id);

num_round = 160;
num_ant = 16;

cb = zeros(num_round,62,num_ant*num_ant);
for i = 1:num_round
    curcb = kron(tx_codebook(1+62*(i-1):62+62*(i-1),:),rx_codebook(i,:));
    cb(i,:,:) = curcb;
end

cb = permute(cb, [2,1,3]);
cb = reshape(cb,[num_round*62,num_ant*num_ant]);

save('./codebook_mat/random_probe_cb_16x16_multires.mat','cb');