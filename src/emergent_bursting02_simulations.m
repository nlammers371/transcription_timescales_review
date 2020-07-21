% conduct small-scale simulations to generate illustrative traces for
% Figure 3
clear
close all
addpath('utilities')

% load numeric results
project = 'n6';
addpath('utilities')

% set paths
DataPath = ['../out/emergent_bursting/' project '/'];

load([DataPath 'bursting_chain_calc_struct.mat'])
load([DataPath 'bursting_step_calc_struct.mat'])

% set stochastic simulation parameters
n_sim = 25; % number of indivudal traces to simulate
t_sim = 60*60; % duration of each simulation (in seconds)
% n_calc_points = size(bursting_chain_calc_struct(1).Q,3);

p = gcp('nocreate');
if isempty(p)
  parpool(24);
end

%%%%%%%%%%%%%%%%%%%%%%%
% independent burst chain  first
%%%%%%%%%%%%%%%%%%%%%%%
ind_index = 1;
ind_sub_index = 101;
r_seed = 123;

bursting_sim_struct = struct;
% set seed for consistency

% record key info

disp('Simulating independent binding condition...')
tic  
bursting_temp = stochastic_simulation_wrapper(bursting_chain_calc_struct(ind_index),ind_sub_index,n_sim,t_sim,r_seed);
fnames = fieldnames(bursting_temp);
for f = 1:length(fnames)
  bursting_sim_struct(1).(fnames{f}) = bursting_temp.(fnames{f});
end
toc

%%%%%%%%%%%%%%%%%%%%%%%
% cooperative binding
%%%%%%%%%%%%%%%%%%%%%%%

% now cooperative chain
simIndicesCoop = 176:4:201;

% set seed for consistency
r_seed_vec = round(rand(1,length(bursting_step_calc_struct)+2)*1000);
disp('Simulating cooperative binding condition...')
tic
iter = 1;
for i = [2 3]
  bursting_temp = stochastic_simulation_wrapper(bursting_chain_calc_struct(i),simIndicesCoop,n_sim,t_sim,r_seed_vec(iter));
  fnames = fieldnames(bursting_temp);
  for f = 1:length(fnames)
    bursting_sim_struct(i).(fnames{f}) = bursting_temp.(fnames{f});
  end
  iter = iter + 1;
end
toc

% now iterate through compound chains for rate limiting step
rng(122);
r_seed_rl_vec = round(rand(1,length(bursting_step_calc_struct))*1000);
simIndicesRL = simIndicesCoop;%[simIndicesCoop(end-2) simIndicesCoop(end)];
offset = length(bursting_sim_struct);

disp('Simulating rate-limiting step codition...')
tic
for i = 2%1:length(bursting_step_calc_struct)
  bursting_temp = stochastic_simulation_wrapper(bursting_step_calc_struct(i),simIndicesRL,n_sim,t_sim,r_seed_rl_vec(iter));
  fnames = fieldnames(bursting_temp);
  for f = 1:length(fnames)
    bursting_sim_struct(iter+1).(fnames{f}) = bursting_temp.(fnames{f});
  end
  iter = iter + 1;
end
toc

save([DataPath, 'bursting_sim_struct.mat'],'bursting_sim_struct');