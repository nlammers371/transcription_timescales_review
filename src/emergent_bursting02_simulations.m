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
n_sim = 50; % number of indivudal traces to simulate
t_sim = 60*60; % duration of each simulation (in seconds)
n_calc_points = size(bursting_chain_calc_struct(1).Q,3);

% independent burst chain  first
bursting_sim_struct = struct;
simIndicesIndependent = [1 101 201];

% set seed for consistency
rng(123)
% record key info
bursting_sim_struct(1).name = bursting_chain_calc_struct(1).name;
bursting_sim_struct(1).E = bursting_chain_calc_struct(1).E;

disp('Simulating independent binding condition...')
tic
for p = 1:length(simIndicesIndependent)
  % record network characteristics
  bursting_sim_struct(1).Q(:,:,p) = bursting_chain_calc_struct(1).Q(:,:,simIndicesIndependent(p));
  bursting_sim_struct(1).SS(:,p) = bursting_chain_calc_struct(1).SS(:,simIndicesIndependent(p));   
  bursting_sim_struct(1).SSFull(:,p) = bursting_sim_struct(1).SS(:,p);
  % call sim function
  
  [~, bursting_sim_struct(1).sim_emission_cell(p,:),bursting_sim_struct(1).sim_time_cell(p,:)] = ...
        stochastic_sim_fun(bursting_sim_struct(1).Q(:,:,p), bursting_sim_struct(1).SS(:,p),bursting_sim_struct(1).E,n_sim,t_sim);  
end
toc

% now cooperative chain
simIndicesCoop = 180:6:201;

% set seed for consistency
rng(231)
disp('Simulating cooperative binding condition...')
tic
for i = [2 3]
  % record key info
  bursting_sim_struct(i).name = bursting_chain_calc_struct(i).name;
  bursting_sim_struct(i).E = bursting_chain_calc_struct(i).E;
  for p = 1:length(simIndicesCoop)
    % record network characteristics
    bursting_sim_struct(i).Q(:,:,p) = bursting_chain_calc_struct(i).Q(:,:,simIndicesCoop(p));
    bursting_sim_struct(i).SS(:,p) = bursting_chain_calc_struct(i).SS(:,simIndicesCoop(p));   
    bursting_sim_struct(i).SSFull(:,p) = bursting_sim_struct(i).SS(:,p);
    % call sim function
   
    [~, bursting_sim_struct(i).sim_emission_cell(p,:),bursting_sim_struct(i).sim_time_cell(p,:)] = ...
          stochastic_sim_fun(bursting_sim_struct(i).Q(:,:,p), bursting_sim_struct(i).SS(:,p),bursting_sim_struct(i).E,n_sim,t_sim);
   
  end
end
toc
% now iterate through compound chains for rate limiting step
offset = length(bursting_sim_struct);
rng(234)
disp('Simulating rate-limiting step codition...')
tic
for i = 1:length(bursting_step_calc_struct)
  % record key info
  bursting_sim_struct(i+offset).name = bursting_step_calc_struct(i).name;
  bursting_sim_struct(i+offset).E = bursting_step_calc_struct(i).E;
  for p = 1:length(simIndicesCoop)
    % record network characteristics
    bursting_sim_struct(i+offset).Q(:,:,p) = bursting_step_calc_struct(i).Q(:,:,simIndicesCoop(p));
    bursting_sim_struct(i+offset).SS(:,p) = bursting_step_calc_struct(i).SS(:,simIndicesCoop(p));   
    bursting_sim_struct(i+offset).SSFull(:,p) = bursting_step_calc_struct(i).SSFull(:,simIndicesCoop(p));   
    % call sim function
    [bursting_sim_struct(i+offset).sim_state_cell(p,:), bursting_sim_struct(i+offset).sim_emission_cell(p,:),bursting_sim_struct(i+offset).sim_time_cell(p,:)] = ...
          stochastic_sim_fun(bursting_sim_struct(i+offset).Q(:,:,p), bursting_sim_struct(i+offset).SSFull(:,p),bursting_sim_struct(i+offset).E,n_sim,t_sim);
  end
end
toc

save([DataPath, 'bursting_sim_struct.mat'],'bursting_sim_struct');