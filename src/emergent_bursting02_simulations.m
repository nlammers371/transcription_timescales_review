% script to simulate burst time series with parameters from analytic
% calculations
clear
close all
addpath('utilities')

% load numeric results
project = 'n6';
addpath('utilities')

% set paths
DataPath = ['../out/emergent_bursting/' project '/'];

% load data
load([DataPath 'bursting_chain_calc_struct.mat'])
load([DataPath 'bursting_step_calc_struct.mat'])

% set stochastic simulation parameters
n_sim = 100; % number of indivudal traces to simulate
t_sim = 60*60; % duration of each simulation (in seconds)
n_calc_points = size(bursting_chain_calc_struct(1).Q,3);

% independent burst chain  first
bursting_sim_struct = struct;
simIndicesIndependent = [1 51 101 151 201];

% set seed for consistency
rng(123)
% record key info
bursting_sim_struct(1).name = bursting_chain_calc_struct(1).name;
bursting_sim_struct(1).E = bursting_chain_calc_struct(1).E;
for p = 1:length(simIndicesIndependent)
  % record network characteristics
  bursting_sim_struct(1).Q(:,:,p) = bursting_chain_calc_struct(1).Q(:,:,simIndicesIndependent(p));
  bursting_sim_struct(1).SS(:,p) = bursting_chain_calc_struct(1).SS(:,simIndicesIndependent(p));   
  % call sim function
  tic
  [bursting_sim_struct(1).sim_state_cell(p,:), bursting_sim_struct(1).sim_emission_cell(p,:),bursting_sim_struct(1).sim_time_cell(p,:)] = ...
        stochastic_sim_fun(bursting_sim_struct(1).Q(:,:,p), bursting_sim_struct(1).SS(:,p),bursting_sim_struct(1).E,n_sim,t_sim);
  toc
end

%% now cooperative chain
simIndicesCoop = 180:3:201;

% set seed for consistency
rng(123)

for i = [2 3]
  % record key info
  bursting_sim_struct(i).name = bursting_chain_calc_struct(1).name;
  bursting_sim_struct(i).E = bursting_chain_calc_struct(1).E;
  for p = 1:length(simIndicesCoop)
    % record network characteristics
    bursting_sim_struct(i).Q(:,:,p) = bursting_chain_calc_struct(i).Q(:,:,simIndicesCoop(p));
    bursting_sim_struct(i).SS(:,p) = bursting_chain_calc_struct(i).SS(:,simIndicesCoop(p));   
    % call sim function
    tic
    [bursting_sim_struct(i).sim_state_cell(p,:), bursting_sim_struct(i).sim_emission_cell(p,:),bursting_sim_struct(i).sim_time_cell(p,:)] = ...
          stochastic_sim_fun(bursting_sim_struct(i).Q(:,:,p), bursting_sim_struct(i).SS(:,p),bursting_sim_struct(i).E,n_sim,t_sim);
    toc
  end
end
%%
save([DataPath, 'bursting_sim_struct.mat'],'bursting_sim_struct');