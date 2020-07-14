% script to simulate burst time series with parameters from analytic
% calculations
function bursting_sim_struct = stochastic_simulation_wrapper(bursting_calc_struct,simIndex,simSubIndices,n_sim,t_sim,r_seed)

addpath('utilities')

% load numeric results
% project = 'n6';
addpath('utilities')

% set paths
% DataPath = ['../out/emergent_bursting/' project '/'];


% set stochastic simulation parameters
% n_sim = 100; % number of indivudal traces to simulate
% t_sim = 60*60; % duration of each simulation (in seconds)
% n_calc_points = size(bursting_calc_struct(1).Q,3);

% independent burst chain  first
bursting_sim_struct = struct;
% set seed for consistency
rng(r_seed)

% record key info
bursting_sim_struct.name = bursting_calc_struct(simIndex).name;
bursting_sim_struct.E = bursting_calc_struct(simIndex).E;
for p = 1:length(simSubIndices)
  % record network characteristics
  bursting_sim_struct.Q(:,:,p) = bursting_calc_struct(simIndex).Q(:,:,simSubIndices(p));
  bursting_sim_struct.SS(:,p) = bursting_calc_struct(simIndex).SS(:,simSubIndices(p));   
  if isfield(bursting_calc_struct,'SSFull')
    bursting_sim_struct.SSFull(:,p) = bursting_calc_struct(simIndex).SSFull(:,p);
  else
    bursting_sim_struct.SSFull(:,p) = bursting_sim_struct.SS(:,p);
  end
  
  % call sim function
  [bursting_sim_struct.sim_emission_cell(p,:), bursting_sim_struct.sim_emission_cell(p,:),bursting_sim_struct.sim_time_cell(p,:)] = ...
        stochastic_sim_fun(bursting_sim_struct.Q(:,:,p), bursting_sim_struct.SSFull(:,p),bursting_sim_struct.E,n_sim,t_sim);

end