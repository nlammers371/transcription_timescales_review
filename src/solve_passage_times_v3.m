% script to solve for effective off/on rates in linear chain binding model
clear
close all
addpath('utilities')

% specify model characteristics
n_bound_vec = 0:11;
n_states = numel(n_bound_vec);

% ON and OFF blocks
on_states = 6:11;
off_states = n_bound_vec(~ismember(n_bound_vec,on_states));
midpoint = 5.5;

% specify binding dynamics
syms kon koff amp n_hill
nonlin_rate = 'kon';
on_nonlin_flag = strcmp(nonlin_rate,'kon');
off_nonlin_flag = strcmp(nonlin_rate,'koff');

% koff_val (from Mir et al 2018)
koff_val = 1/3;

% define null rates without nonlinearity
koff_vec = koff*n_bound_vec; % assume Bcd unbinds independently
kon_vec = kon*(n_states-n_bound_vec-1); % assume Bcd unbinds independently

if on_nonlin_flag    
    kon_vec = kon_vec.*(1+amp .* n_bound_vec.^n_hill ./ (midpoint.^n_hill + n_bound_vec.^n_hill));
end
if off_nonlin_flag   
    koff_vec = koff_vec.*(1+amp*(1 - n_bound_vec.^n_hill ./ (midpoint.^n_hill + n_bound_vec.^n_hill)));
end

% state dwell time vec
tau_vec = 1./(kon_vec + koff_vec);

% specify simulation type
project = ['n' num2str(max(n_bound_vec)) '_' nonlin_rate];
% set figure path
FigPath = ['../fig/emergent_bursting/' project '/'];
mkdir(FigPath)
DataPath = ['../out/emergent_bursting/' project '/'];
mkdir(DataPath)

%% generate expressions for first passage times from ON to OFF blocks
% OFF to ON
off_eq_list = {};
off_sym_list = [];
off_str_list = {};

% generate list of symbols to use
for i = off_states
    % initialize symbol
    sym_str = ['ET' num2str(i) num2str(on_states(1))];
    off_str_list{i+1} = sym_str;
    eval(['syms ' sym_str])
    off_sym_list = [off_sym_list eval(sym_str)];
end

for i = (off_states)+1
    % generate equation
    next_on = 0;
    if i < max(off_states)+1
        next_on = off_sym_list(i+1);
    end
    next_off = 0;
    if i > min(off_states)+1
        next_off = off_sym_list(i-1);
    end
    off_eq_list{i} = off_sym_list(i) == tau_vec(i)*(1 + kon_vec(i)*next_on + koff_vec(i)*next_off);
end
% solve 
solOff = solve(off_eq_list,off_sym_list);
eval(['off_dwell = solOff.' off_str_list{end} ';']);
kon_eff_fun = matlabFunction(1/off_dwell);
kon_sym = 1/off_dwell;


% ON to OFF
on_eq_list = {};
on_sym_list = [];
on_str_list = {};
% generate list of symbols to use
for i = on_states
    % initialize symbol
    sym_str = ['ET' num2str(i) num2str(off_states(end))];
    on_str_list = [on_str_list{:} {sym_str}];
    eval(['syms ' sym_str])
    on_sym_list = [on_sym_list eval(sym_str)];
end
%
os = numel(off_states);
for i = 1:numel(on_states)
    % generate equation
    next_on = 0;
    if i < numel(on_states)
        next_on = on_sym_list(i+1);
    end
    next_off = 0;
    if i > 1
        next_off = on_sym_list(i-1);
    end
    on_eq_list = [on_eq_list {on_sym_list(i) == tau_vec(i+os)*(1 + kon_vec(i+os)*next_on + koff_vec(i+os)*next_off)}];
end
% solve 
solOn = solve(on_eq_list,on_sym_list);
eval(['on_dwell = solOn.' on_str_list{1} ';']);
koff_eff_fun = matlabFunction(1/on_dwell);


%% conduct simple parameter sweep of possible solutions'
n_samp = 1e4;
kon_samples = 10.^(rand(1,n_samp)*-3);
n_hill_samples = 10.^(rand(1,n_samp)*2);
amp = 100;

kon_eff_vec = NaN(1,n_samp);
koff_eff_vec = NaN(1,n_samp);
tic
for i = 1:n_samp
    kon_eff_vec(i) = 60*kon_eff_fun(amp,koff_val,kon_samples(i),n_hill_samples(i));%(kon_samples(i),1,b_samples(i)));
    koff_eff_vec(i) = 60*koff_eff_fun(amp,koff_val,kon_samples(i),n_hill_samples(i));%double(koff_eff(kon_samples(i),1,b_samples(i)));
end
toc

slow_ft = kon_eff_vec <= 1 & koff_eff_vec <= 1;
cmap = brewermap([],'Set2');
close all

scatter_fig = figure;
hold on
scatter(kon_eff_vec,koff_eff_vec,20,cmap(3,:))
scatter(kon_eff_vec(slow_ft),koff_eff_vec(slow_ft),20,cmap(2,:))
set(gca,'YScale','log','XScale','log')    
set(gca,'Fontsize',14)
xlabel('effective on rate (events/min)')
ylabel('effective off rate (events/min)')
ylim([1e-4 1e2])
grid on    
saveas(scatter_fig,[FigPath 'param_sweep.png'])    
    
% save analytic results
results_struct = struct;
results_struct.kon_eff_fun = kon_eff_fun;
results_struct.koff_eff_fun = koff_eff_fun;
results_struct.koff_val = koff_val;

results_struct.kon_eff_vec = kon_eff_vec;
results_struct.koff_eff_vec = koff_eff_vec;
results_struct.kon_samples = kon_samples;
results_struct.b_samples = amp_samples;

results_struct.n_bound_vec = n_bound_vec;
results_struct.on_states = on_states;
results_struct.off_states = off_states;
results_struct.kon_vec = kon_vec;
results_struct.koff_vec = koff_vec;
results_struct.tau_vec = tau_vec;

save([DataPath 'results_struct.mat'],'results_struct')