% Run concurrent process parameter estimation

clear persistent
clear

for i =1
[num,txt,raw] = xlsread('params_par.xlsx');
num_in=1;

% Parameter bounds
lb_cd = 7; ub_cd = 13;
lb_bii = 4.5; ub_bii = 6;
lb_bbd = 1.1; ub_bbd = 2.2;
lb_sigii=0.05; ub_sigii= 0.4;
lb_sigcd= 0.1; ub_sigcd = 0.6;
lb_sigbbd= 0.01; ub_sigbbd = 0.5;

lbs = [lb_cd, lb_bii, lb_bbd, lb_sigii, lb_sigcd, lb_sigbbd];
ubs = [ub_cd, ub_bii, ub_bbd, ub_sigii, ub_sigcd, ub_sigbbd];

[td, lb, ld, li, cd, gr, moth_dau] = analyse_cycle_data();
observed_data = [td, lb, ld, li, cd, (ld-lb)./td];

% Optimization options
options = optimoptions('ga', 'Display', 'final', 'UseParallel',true,'MaxGenerations',50,'PlotFcn',"gaplotbestf",'OutputFcn', @gaOutputFcn);

lossWrapper = @(params) concurrent_process_bir_add_init_timer_cond_prob(params, observed_data, num, num_in);

% Run the optimization
[param_opt, loss_opt] = ga(lossWrapper, 6, [], [], [], [], lbs, ubs, [], options);

disp('Optimal Parameters:');
disp(param_opt);
disp('Optimal Loss:');
disp(loss_opt);

end


%%

% Run CH model parameter estimation

clear persistent
clear

for i = 1

[num,txt,raw] = xlsread('params_par.xlsx');
num_in=1;

% Parameter bounds
lb_cd = 11; ub_cd = 15; 
lb_bii= 4.8; ub_bii = 6;
lb_sigii= 0.05; ub_sigii = 0.45;
lb_sigcd= 0.15; ub_sigcd = 0.6;

lbs = [lb_cd, lb_bii, lb_sigii, lb_sigcd];
ubs = [ub_cd, ub_bii, ub_sigii, ub_sigcd];

[td, lb, ld, li, cd, gr, moth_dau] = analyse_cycle_data();
observed_data = [td, lb, ld, li, cd, (ld-lb)./td];

% Optimization options
options = optimoptions('ga', 'Display', 'final', 'UseParallel',true,'MaxGenerations',50,'PlotFcn',"gaplotbestf",'OutputFcn', @gaOutputFcn);

lossWrapper = @(params) adder_per_origin(params, observed_data, num, num_in);

% Run the optimization
[param_opt, loss_opt] = ga(lossWrapper, 4, [], [], [], [], lbs, ubs, [], options);

disp('Optimal Parameters:');
disp(param_opt);
disp('Optimal Loss:');
disp(loss_opt);
end

%%

% Run parallel adder optimization

clear persistent
clear

for i = 1
[num,txt,raw] = xlsread('params_par.xlsx');
num_in=1;

% Parameter bounds
lb_bid = 1.4; ub_bid = 2.4;
lb_bii= 4.6; ub_bii = 6;
lb_sigii= 0.05; ub_sigii = 0.45;
lb_sigid= 0.1; ub_sigid = 0.6;

lbs = [lb_bid, lb_bii, lb_sigii, lb_sigid];
ubs = [ub_bid, ub_bii, ub_sigii, ub_sigid];

[td, lb, ld, li, cd, gr, moth_dau] = analyse_cycle_data();
observed_data = [td, lb, ld, li, cd, (ld-lb)./td];

% Optimization options
options = optimoptions('ga', 'Display', 'final', 'UseParallel',true,'MaxGenerations',50,'PlotFcn',"gaplotbestf",'OutputFcn', @gaOutputFcn);

lossWrapper = @(params) parallel_adder(params, observed_data, num, num_in);

% Run the optimization
[param_opt, loss_opt] = ga(lossWrapper, 4, [], [], [], [], lbs, ubs, [], options);

disp('Optimal Parameters:');
disp(param_opt);
disp('Optimal Loss:');
disp(loss_opt);

end
