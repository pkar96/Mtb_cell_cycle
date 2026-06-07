% calls adder per origin/ parallel adder model /concurrent process 
% and get simulation results stored in variables similar to that in 
% experiments. Call analyse_cycle_data, cv_plots and gen_sim_plots
% to generate plots.

% Steps to follow 1: Uncomment the section for the plots to do- 
% CV or correlation.

clear

% Run simulations - choose model

model_n = 0; % 0 - experimental data, model_n values for the models are specified below

    if model_n == 1
        [tds, vbs, vds, vios, cds,rts, moth_infos] = adder_per_origin([11.8219 5.7676 0.2725 0.31577]); % CH model
    elseif model_n == 2
        [tds, vbs, vds, vios, cds,rts, moth_infos] = parallel_adder([1.7123 5.4513 0.2062 0.3363]); % parallel adder model
    elseif model_n == 9
        [tds, vbs, vds, vios, cds,rts, moth_infos] = concurrent_process_bir_add_init_timer_cond_prob([10.362 5.3573 1.8525 0.20843 0.24793 0.32974]); % concurrent processes
    end

% Run exp- choose media
[td, lb, ld, li, cd, gr, moth_dau] = analyse_cycle_data();
%%
%-----plotting CV values
if model_n==0
    figure
    hold on
    cv_plots(td, lb, ld, li, cd, gr, 0)
else
    cv_plots(tds, vbs, vds, vios, cds, rts, model_n)
end
set(gca, 'FontSize', 24)
set(gcf, 'Position',[276,42,777,602])
%%
%-----plotting correlation values
if model_n==0
    figure
    hold on
    gen_sim_plots(td, lb, ld, li, cd, gr, moth_dau, 0)
else
    gen_sim_plots(tds, vbs, vds, vios, cds, rts, moth_infos, model_n)
end
set(gca, 'FontSize', 24)
set(gcf, 'Position',[276,42,777,602])