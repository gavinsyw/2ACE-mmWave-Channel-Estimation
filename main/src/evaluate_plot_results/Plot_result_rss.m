function Plot_result_rss(Label_Name, Simulation_result)

plot_init;

%% NMSE of Channel Matrix
figure;
%-----NonCoherent/PerfectPhase/NoisePhase-----
[v, Active_ind] = find(Simulation_result.Method.State>0);
for i = Active_ind
    plot(Simulation_result.Range,10*log10(Simulation_result.Mean_Evaluation(1,:,4,i))); 
    hold on;
end
%-----Labels and legends-----
xlabel(Label_Name,'Interpreter','latex');
ylabel('NMSE of Channel Matrix H (dB)','Interpreter','latex');  
legend_name = Simulation_result.Method.Name_of_method;
legend_name = legend_name(Active_ind);        
legend(legend_name); grid on;
if contains(Label_Name,'Searching') 
    set(gca, 'XDir','reverse')
end

end