% script to simulate burst time series with parameters from analytic
% calculations
function bursting_sim_struct = stochastic_simulations_wrapper(DataPath,simIndex,simSubIndices,n_sim,t_sim,r_seed)
% clear
% close all
addpath('utilities')

% load numeric results
% project = 'n6';
addpath('utilities')

% set paths
% DataPath = ['../out/emergent_bursting/' project '/'];

% load data
load([DataPath 'bursting_chain_calc_struct.mat'])
load([DataPath 'bursting_step_calc_struct.mat'])

% set stochastic simulation parameters
% n_sim = 100; % number of indivudal traces to simulate
% t_sim = 60*60; % duration of each simulation (in seconds)
% n_calc_points = size(bursting_chain_calc_struct(1).Q,3);

% independent burst chain  first
bursting_sim_struct = struct;

% now cooperative chain
% simSubIndices = 176:4:201;

% set seed for consistency
rng(r_seed)
disp('Simulating cooperative binding condition...')
tic

% record key info
bursting_sim_struct(simIndex).name = bursting_chain_calc_struct(simIndex).name;
bursting_sim_struct(simIndex).E = bursting_chain_calc_struct(simIndex).E;
for p = 1:length(simSubIndices)
  % record network characteristics
  bursting_sim_struct(simIndex).Q(:,:,p) = bursting_chain_calc_struct(simIndex).Q(:,:,simSubIndices(p));
  bursting_sim_struct(simIndex).SS(:,p) = bursting_chain_calc_struct(simIndex).SS(:,simSubIndices(p));   
  bursting_sim_struct(simIndex).SSFull(:,p) = bursting_sim_struct(simIndex).SS(:,p);
  % call sim function

  [bursting_sim_struct(simIndex).sim_emission_cell(p,:), bursting_sim_struct(simIndex).sim_emission_cell(p,:),bursting_sim_struct(simIndex).sim_time_cell(p,:)] = ...
        stochastic_sim_fun(bursting_sim_struct(simIndex).Q(:,:,p), bursting_sim_struct(simIndex).SS(:,p),bursting_sim_struct(simIndex).E,n_sim,t_sim);

end
toc

% save([DataPath, 'bursting_sim_struct.mat'],'bursting_sim_struct');