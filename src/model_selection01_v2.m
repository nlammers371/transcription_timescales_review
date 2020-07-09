% script to examine waiting time distributions for bottleneck and emergent
% bursting simulations
clear
close all
addpath('utilities')

% load numeric results
n_bcd_sites = 6;
project = ['n' num2str(n_bcd_sites)];
addpath('utilities')

% set paths
FigurePath = ['../fig/emergent_bursting/' project '/'];
mkdir(FigurePath)
DataPath = ['../out/emergent_bursting/' project '/'];

% load data
load([DataPath 'bursting_sim_struct.mat'])
load([DataPath 'bursting_chain_calc_struct.mat'])

% set basic parameters

emergent_indices = 2:3;
rate_lim_indices = 4:length(bursting_sim_struct);
sim_indices = [emergent_indices rate_lim_indices];
dT = 1; % time res for interpolated data
time_vector = 0:dT:60*60;
n_bound_vec = 0:6;
n_bs = 6;
sim_index = 8;
%% (1) fit 2 state HMM to simulated data
n_sim = size(bursting_sim_struct(1).sim_time_cell,2);

% we first need to interpolate the simulation results
trace_array = NaN(n_sim, length(time_vector),length([emergent_indices rate_lim_indices]));

for i = 1:size(trace_array,3)
  for n = 1:n_sim
    % emergent traces first 
    trace = bursting_sim_struct(sim_indices(i)).sim_emission_cell{sim_index,n}*n_bs + 1;
    time = bursting_sim_struct(sim_indices(i)).sim_time_cell{sim_index,n};
    trace_array(n,:,i) = interp1(time,trace,time_vector,'previous','extrap');
  end
end

%% now plug traces into HMM

% initial guess for transition prob matrix
A_guess = ones(2);
A_guess(eye(2)==1) = 50;
A_guess = A_guess./sum(A_guess);

% guess for emission probabilities
E_guess = ones(2,length(n_bound_vec));
E_guess(1,1) = 10;
E_guess(2,end) = 10;
E_guess = E_guess ./ sum(E_guess,2);

% initialize arrays 
A_array = NaN(2,2,size(trace_array,3));
E_array = NaN(2,n_bs+1,size(trace_array,3));

% estimate HMM for the two models 
disp('estimating HMM models...')
for i = 1:size(trace_array,3)
  tic
  [A_array(:,:,i),E_array(:,:,i)] = hmmtrain(trace_array(:,:,i),A_guess,E_guess);
  toc
end

% estimate viterbi paths
viterbi_traces = NaN(size(trace_array));

disp('performing Viterbi fits...')
for i = 1:size(trace_array,3)
  tic
  for n = 1:n_sim
    viterbi_traces(n,:,i) = hmmviterbi(trace_array(n,:,i),A_array(:,:,i),E_array(:,:,i));
  end
  toc
end


%% (2) Now extract waiting times for each
wt_off_cell = cell(1,length(sim_indices));

for i = 1:length(sim_indices)
  wt_off_vec = [];

  for n = 1:n_sim
    % emergent
    ev_trace = viterbi_traces(n,:,i);
    ev_diff = [0 diff(ev_trace)];
    d_points = find(ev_diff);
    dt_vec = diff(d_points)*dT;
    id_vec = diff(ev_diff(d_points));
    
    wt_off_vec = [wt_off_vec dt_vec(id_vec>0)];   

  end
  wt_off_cell{i} = wt_off_vec;
end
%%
e_vec = [0 0 1];
Q = ones(3);
Q(eye(3)==1) = 0;
Q(3,1) = 0;
Q(1,2) = 0;
Q(2,1) = 1;
Q(2,3) = 0;
Q(eye(3)==1) = -sum(Q);
state_options = 1:3;
state_vec = [1];
time_vec = [0];
T = 1e5;
tc = 0;
while tc < T
  state_curr = state_vec(end);
  dt = exprnd(-1/Q(state_curr,state_curr));
  tc = tc + dt;
  
  w_vec = Q(:,state_curr);
  w_vec(state_curr) = 0;
  next_state = randsample(state_options,1,true,w_vec);
  
  if tc < T
    state_vec(end+1) = next_state;
    time_vec(end+1) = tc;
  end
end
%%
e_vec = interp1(time_vec,1*(state_vec == 3),linspace(0,max(time_vec),10*T),'prev');
e_diff = [0 diff(e_vec)];
d_points = find(e_diff);
dt_vec = diff(d_points)*dT;
id_vec = diff(e_diff(d_points));