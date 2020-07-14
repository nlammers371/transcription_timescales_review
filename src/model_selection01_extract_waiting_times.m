% script to generate data for wiating time distributions of different
% models
clear
close all
addpath('utilities')

% load numeric results
n_bcd_sites = 6;
project = ['n' num2str(n_bcd_sites)];
addpath('utilities')

% set paths
ReadPath = ['../out/emergent_bursting/' project '/'];
WritePath = ['../out/waiting_time_distributions/' project '/'];
mkdir(WritePath);

% load data
load([ReadPath 'bursting_sim_struct.mat'])
load([ReadPath 'bursting_chain_calc_struct.mat'])

% set basic parameters

emergent_indices = 2:3;
rate_lim_indices = 4:length(bursting_sim_struct);
sim_indices = [emergent_indices rate_lim_indices];

dT = 1; % time res for interpolated data in seconds

time_vector = 0:dT:60*60;
n_bound_vec = 0:6;
n_bs = 6;

waiting_time_struct = struct;
iter = 1;
for s = sim_indices
  
  sub_index_vec = 1:size(bursting_sim_struct(s).SS,2);

  % (1) fit 2 state HMM to simulated data
  n_sim = size(bursting_sim_struct(1).sim_time_cell,2); % number of independent simulations

  % we first need to interpolate the simulation results
  trace_array = NaN(n_sim, length(time_vector),length(sub_index_vec));

  for i = 1:size(trace_array,3)
    for n = 1:n_sim
      % emergent traces first 
      trace = single(bursting_sim_struct(s).sim_emission_cell{i,n} + 1);
      time = bursting_sim_struct(s).sim_time_cell{i,n};
      trace_array(n,:,i) = interp1(time,trace,time_vector,'previous','extrap');
    end
  end

  % now plug traces into HMM

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
  disp(['estimating HMM models (' num2str(iter) '/' num2str(length(sim_indices)) ')...'])
  rng(346);
  tic
  parfor i = 1:size(trace_array,3)    
    [A_array(:,:,i),E_array(:,:,i)] = hmmtrain(trace_array(:,:,i),A_guess,E_guess);    
  end
  toc

  % estimate viterbi paths
  viterbi_traces = NaN(size(trace_array));

%   disp('performing Viterbi fits...')
%   tic
  for i = 1:size(trace_array,3)    
    for n = 1:n_sim
      viterbi_traces(n,:,i) = hmmviterbi(trace_array(n,:,i),A_array(:,:,i),E_array(:,:,i))-1;
    end
  end
%   toc


  % (2) Now extract waiting times for each
  wt_off_cell = cell(1,length(sim_indices));

  for i = 1:length(sub_index_vec)
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

  % store results in a data structure  
  waiting_time_struct(iter).A_array = A_array;
  waiting_time_struct(iter).E_array = E_array;
  waiting_time_struct(iter).trace_array = trace_array;
  waiting_time_struct(iter).time_vector = time_vector;
  waiting_time_struct(iter).viterbi_traces = viterbi_traces;
  waiting_time_struct(iter).off_waiting_times = wt_off_cell;
  waiting_time_struct(iter).name = bursting_sim_struct(s).name;
  iter = iter + 1;
end

% save
save([WritePath 'waiting_time_struct.mat'],'waiting_time_struct')