clear
clc

load('../result/LoS_2023-08-18.mat');

tx_ant_num = 16;
rx_ant_num = 16;
Mt = round(linspace(2,sqrt(4*tx_ant_num*rx_ant_num),8));
Mr = Mt;

data = zeros(8,10);
for i = 1:8
    for j = 1:10
        if j == 1
            data(i,j) = max(squeeze(rss_bf(:, i,j, 1:50)),[],'all');
        elseif j == 2
            M_cur = Mt(i);
            data(i,j) = beam_sweeping(rss_sweeping_phi,M_cur,32, 10000);
        elseif j == 3
            M_cur = Mt(i);
            data(i,j) = beam_sweeping(rss_sweeping_thetaNphi,M_cur, 36,10000);
        elseif j == 4
            temp = sort(max(squeeze(rss_bf(:, i,j-2,51:52)),[],2),'descend');
            data(i,j) = max(temp,[],'all');
            %data(i,j) = mean(temp(1:5));
        elseif j == 5
            temp = sort(max(squeeze(rss_bf(:, i,j-2,53:54)),[],2),'descend');
            data(i,j) = max(temp,[],'all');
            %data(i,j) = mean(temp(1:5));
        elseif j == 6
            temp = sort(max(squeeze(rss_bf(:, i,j-2,55:56)),[],2),'descend');
            data(i,j) = max(temp,[],'all');
            %data(i,j) = mean(temp(1:5));
        elseif j == 7
            temp = sort(max(squeeze(rss_bf(:, i,j-2,57:58)),[],2),'descend');
            data(i,j) = max(temp,[],'all');
            %data(i,j) = mean(temp(1:5));
        elseif j == 8
            temp = sort(max(squeeze(rss_bf(:, i,j-2,59:60)),[],2),'descend');
            data(i,j) = max(temp,[],'all'); 
            %data(i,j) =  mean(temp(1:5));
        elseif j == 9
            temp = sort(max(squeeze(rss_bf(:, i,j-2,61:62)),[],2),'descend');
            data(i,j) = max(temp,[],'all');
            %data(i,j) = mean(temp(1:5));
        elseif j == 10
            temp = sort(max(squeeze(rss_bf(:, i,j-2,63:64)),[],2),'descend');
            data(i,j) = max(temp,[],'all');
            %data(i,j) = mean(temp(1:5));
        end
    end
end

data(:,1) = mean(data(:,1));
data(:,10) = mean(data(:,10));

x = [4,36,121,225,361,529,784,1024];
createfigure(x,data,[-63,-53])

function rss = beam_sweeping(rss_in, tx_beam, total_beam, num_run)
    rss = zeros(1,num_run);
    for i = 1:num_run
        idx = randperm(total_beam,tx_beam);
        rss(i) = max(rss_in(idx,idx),[],"all");
    end
    rss = mean(rss);
end