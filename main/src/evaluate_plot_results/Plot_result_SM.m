%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2019 Yi Zhang and The University of Texas at Austin 
%  
% Permission is hereby granted, free of charge, to any person obtaining a 
% copy of this software and associated documentation files (the
% "Software"), to deal in the Software without restriction, including 
% without limitation the rights to use, copy, modify, merge, publish, 
% distribute, sublicense, and/or sell copies of the Software, and to permit
% persons to whom the Software is furnished to do so, subject to the 
% following conditions:
% 
% The above copyright notice and this permission notice shall be included
% in all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
% NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
% DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
% OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
% USE OR OTHER DEALINGS IN THE SOFTWARE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If you use this code or any (modified) part of it in any publication,
% please cite:
%
% Yi Zhang, Kartik Patel, Sanjay Shakkottai, and Robert W. Heath Jr.. 2019. 
% Side-information-aided Non-coherent Beam Alignment Design for Millimeter 
% Wave Systems. In MobiHoc '19: The Twentieth ACM International Symposium 
% on Mobile Ad Hoc Networking and Computing, July 02-05, 2019, Catania, 
% Italy. ACM, New York, NY, USA, 10 pages.
%
% Author: Yi Zhang
% Contact email: yi.zhang.cn@utexas.edu 
% Last modified: Apr. 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function description:
% This function helps to plot the simulation results of the proposed
% algorithm and the related benchmarking algorithms.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% Label_Name: X-axis label (character array).
% Simulation_result: a structure array that groups all related simulation.
% result obtained. Please refer to main scripts in folder main_programs for
% its fields.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Plot_result_SM(Label_Name, Simulation_result)
    plot_init;
    max_ind = numel(Simulation_result.Beampattern_Modes);
    
%     %% Success rate of recovery
%     if Simulation_result.L == 1
%         Num_Quantization_Error = Simulation_result.Num_Quantization_Error;
%         Max_Error = 20;
%         figure
%         j = 1;
%         for i = 1:max_ind
%             Num_Quantization_Error_i = Num_Quantization_Error(:,:,1,i); 
%             Dist_i = hist(Num_Quantization_Error_i,0:Max_Error)/length(Num_Quantization_Error(:,1,1,1));
%             subplot(round(sqrt(max_ind)),ceil(sqrt(max_ind)),j);
%             j = j+1;
%             %colormap hsv;
%             bar3(Dist_i);
%             set(gca,'XTickLabel',Simulation_result.Range);
%             set(gca,'YTickLabel',0:Max_Error)
%             xlabel(Label_Name);
%             ylabel('Number of Quantization Error');
%             zlabel('Percentage');
%             title('Distribution of Recovery In Grid Level '+Simulation_result.Beampattern_Modes(i));
%             view(-50,30);
%             camlight left;
%         end    
%     end
    
    
    %% Error of angle in degree
    figure;
    %-----NonCoherent/PerfectPhase/NoisePhase-----
    for i = 1:max_ind
        plot(Simulation_result.Range,Simulation_result.Mean_Evaluation(1,:,2,i)); 
        hold on;
    end
%     %----Beam sweep----
%     plot(Simulation_result.Range,Simulation_result.Mean_Evaluation(1,:,12,i));
    % label
    xticks(Simulation_result.Range)
    xlabel(Label_Name,'Interpreter','latex');
    ylabel('Mean Angle Estimation Error (MAEE) in Degree','Interpreter','latex');        
    legend_name = Simulation_result.Beampattern_Modes;
%     legend_name = [legend_name, "Beam Sweeping"];
    legend(legend_name,'Interpreter','none'); grid on;
    if contains(Label_Name,'Searching') 
        set(gca, 'XDir','reverse')
    end

    
    %% NMSE of Array Response Vector
    figure;
    %-----NonCoherent/PerfectPhase/NoisePhase-----
    for i = 1:max_ind
        plot(Simulation_result.Range,10*log10(Simulation_result.Mean_Evaluation(1,:,3,i))); 
        hold on;
    end
    %-----Labels and legends-----
    xlabel(Label_Name,'Interpreter','latex');
    ylabel('NMSE of Array Response Vector (dB)','Interpreter','latex'); 
    legend_name = Simulation_result.Beampattern_Modes;
%     legend_name = legend_name(max_ind);        
    legend(legend_name,'Interpreter','none'); grid on;
    if contains(Label_Name,'Searching') 
        set(gca, 'XDir','reverse')
    end


    %% NMSE of Channel Matrix
    figure;
    %-----NonCoherent/PerfectPhase/NoisePhase-----
    for i = 1:max_ind
        plot(Simulation_result.Range,10*log10(Simulation_result.Mean_Evaluation(1,:,4,i))); 
        hold on;
    end
    %-----Labels and legends-----
    xlabel(Label_Name,'Interpreter','latex');
    ylabel('NMSE of Channel Matrix H (dB)','Interpreter','latex');  
    legend_name = Simulation_result.Beampattern_Modes;
%     legend_name = legend_name(max_ind);        
    legend(legend_name,'Interpreter','none'); grid on;
    if contains(Label_Name,'Searching') 
        set(gca, 'XDir','reverse')
    end

    
%     %% Spectrum Efficiency
%     figure;
%     %-----Perfect CSI-----
%     plot(Simulation_result.Range,Simulation_result.Mean_Evaluation(1,:,6,i));hold on;
%     %-----NonCoherent/PerfectPhase/NoisePhase-----
%     for i = 1:max_ind
%         plot(Simulation_result.Range,Simulation_result.Mean_Evaluation(1,:,7,i));
%         hold on;
%     end
%     %-----Beam Sweep-----
%     plot(Simulation_result.Range,Simulation_result.Mean_Evaluation(1,:,10,i));hold on;
%     %-----Labels and legends-----
%     xlabel(Label_Name,'Interpreter','latex');
%     xticks(Simulation_result.Range)
%     ylabel('Spectrum Efficiency (Bits/s/hz)','Interpreter','latex');
%     legend_name = Simulation_result.Beampattern_Modes;
% %     legend_name = legend_name(max_ind);
%     legend_name = ["With Perfect CSI",legend_name,"Beam Sweeping"];            
%     legend(legend_name,'Interpreter','none'); grid on;
%     if contains(Label_Name,'Searching') 
%         set(gca, 'XDir','reverse')
%     end 
end

