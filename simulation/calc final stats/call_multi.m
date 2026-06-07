% calls adder per origin/ parallel adder model /concurrent process 
% and gets AIC values

clear

% observed data
[td, lb, ld, li, cd, gr, moth_dau] = analyse_cycle_data();
observed_data = [td, lb, ld, li, cd, (ld-lb)./td];

% no of sims
nsims = 100;

% net output 
aic_arrs = zeros(nsims, 3);

for nruns = 1:nsims

% output
aic_arr = [];

% Run simulations - choose model
for model_n = [1:2, 9]
    if model_n == 1
        [tds, vbs, vds, vios, cds,rts, moth_infos] = adder_per_origin([11.8219 5.7676 0.2725 0.31577]);
        k=4;
    elseif model_n == 2
        [tds, vbs, vds, vios, cds,rts, moth_infos] = parallel_adder([1.7123 5.4513 0.2062 0.3363]);
        k=4;
    elseif model_n == 9
        [tds, vbs, vds, vios, cds,rts, moth_infos] = concurrent_process_bir_add_init_timer_cond_prob([10.362 5.3573 1.8525 0.20843 0.24793 0.32974]);
        k=6;
    end

    simulated_data = [tds(tds~=0), vbs(tds~=0), vds(tds~=0), vios(tds~=0), cds(tds~=0) rts(tds~=0)];
    aic = aic_data(simulated_data,observed_data,k);

    aic_arr = [aic_arr aic];

end

aic_arrs(nruns,1:length(aic_arr)) = aic_arr;
end

% Display resultsof AIC and Delta AIC
maic = mean(aic_arrs)
maic_d = min(maic(maic~=0))- maic
% aic_arrs_6 = aic_arrs(1:nsims, 3)-aic_arrs(1:nsims, 1:length(aic_arr));

% aic_diff = (min(aic_arr)-aic_arr);
% aic_e = exp(aic_diff/2);
% prob = 1/sum(aic_e);
