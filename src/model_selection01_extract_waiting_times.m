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
load([ReadPath 'bursting_step_calc_struct.mat'])
load([ReadPath 'bursting_chain_calc_struct.mat'])

% simulation parameters
n_sim = 150; % number of indivudal traces to simulate
t_sim = 10*60*60; % duration of each simulation (in seconds)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run stochastic simulation script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = gcp('nocreate');
if isempty(p)
  parpool(24);
end
%%%%%%%%%%%%%%%%%%%%%%%
% cooperative binding
%%%%%%%%%%%%%%%%%%%%%%%

% now cooperative chain
simIndices = 176:4:201;

% set seed for consistency
r_seed_vec = round(rand(1,2+length(bursting_step_calc_struct))*1000);

disp('Simulating cooperative binding condition...')
tic
iter = 1;
for i = [2 3]
  bursting_temp = stochastic_simulation_wrapper(bursting_chain_calc_struct(i),simIndices,n_sim,t_sim,r_seed_vec(iter));
  fnames = fieldnames(bursting_temp);
  for f = 1:length(fnames)
    bursting_sim_struct(iter).(fnames{f}) = bursting_temp.(fnames{f});
  end
  iter = iter + 1;
end
toc

% now iterate through compound chains for rate limiting step
rng(122);

% simIndices = simIndicesCoop;%[simIndicesCoop(end-2) simIndicesCoop(end)];
offset = length(bursting_sim_struct);

disp('Simulating rate-limiting step codition...')
tic
for i = 1:length(bursting_step_calc_struct)
  bursting_temp = stochastic_simulation_wrapper(bursting_step_calc_struct(i),simIndices,n_sim,t_sim,r_seed_vec(i));
  fnames = fieldnames(bursting_temp);
  for f = 1:length(fnames)
    bursting_sim_struct(i+offset).(fnames{f}) = bursting_temp.(fnames{f});
  end
  iter = iter + 1;
end
toc

%%
% set basic parameters
n_bootstraps = 25;
emergent_indices = 1:2;
rate_lim_indices = 3:length(bursting_sim_struct);
sim_indices = [emergent_indices rate_lim_indices];

dT = 5; % time res for interpolated data in seconds (effectively our experimental time res)

time_vector = 0:dT:60*60;
n_bound_vec = 0:6;
n_bs = 6;

waiting_time_struct = struct;
iter = 1;
for s = sim_indices
  
  sub_index_vec = 1:size(bursting_sim_struct(s).SS,2);

  % (1) fit 2 state HMM to simulated data
  n_sim = size(bursting_sim_struct(1).sim_time_cell,2); % number of independent simulations
  
  
  %%% two approaches to extracting waiting times:
  %   (1) HMM: use 2 state HMM to extract state fluctuations. This in
  %       effect mimics what we'd have to do with real data
  
  %   (2) Idealized: extract true times it takes to go from 0->6. 
  
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
  boot_indices = randsample(1:n_sim,n_bootstraps,true); % use subset for speed
  parfor i = 1:size(trace_array,3)    
    [A_array(:,:,i),E_array(:,:,i)] = hmmtrain(trace_array(boot_indices,:,i),A_guess,E_guess);    
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
  wt_off_cell_hmm = cell(1,length(sim_indices));
  wt_off_cell_ideal = cell(1,length(sim_indices));

  for i = 1:length(sub_index_vec)
    wt_off_vec_hmm = [];
    wt_off_vec_ideal = [];
    for n = 1:n_sim
      % extract wiating times using HMM method
      v_trace = viterbi_traces(n,:,i);
      v_diff = [0 diff(v_trace)];
      d_points = find(v_diff);
      dt_vec = diff(d_points)*dT;
      id_vec = diff(v_diff(d_points));

      wt_off_vec_hmm = [wt_off_vec_hmm dt_vec(id_vec>0)];   
      
      % extract idealized waiting times
      raw_trace = bursting_sim_struct(s).sim_emission_cell{i,n};
      raw_times = bursting_sim_struct(s).sim_time_cell{i,n};
      % select for high and low periods
      hl_indices = find(raw_trace==0|raw_trace==n_bcd_sites);     
      condensed_trace = raw_trace(hl_indices);
      % find rise and fall events
      diff_trace = diff([0 condensed_trace]);
      drop_indices = find(diff_trace<0);
      rise_indices = find(diff_trace>0);
      rise_indices = rise_indices(rise_indices>drop_indices(1));
      drop_indices = drop_indices(drop_indices<rise_indices(end));
      % convert back to real time space
      drop_times = raw_times(hl_indices(drop_indices));
      rise_times = raw_times(hl_indices(rise_indices));
      % record
      wt_off_vec_ideal = [wt_off_vec_ideal rise_times-drop_times];
    end
    wt_off_cell_hmm{i} = wt_off_vec_hmm;
    wt_off_cell_ideal{i} = wt_off_vec_ideal;
  end

  % store results in a data structure  
  waiting_time_struct(iter).A_array = A_array;
  waiting_time_struct(iter).E_array = E_array;
  waiting_time_struct(iter).trace_array = trace_array;
  waiting_time_struct(iter).time_vector = time_vector;
  waiting_time_struct(iter).viterbi_traces = viterbi_traces;
  waiting_time_struct(iter).off_waiting_times_hmm = wt_off_cell_hmm;
  waiting_time_struct(iter).off_waiting_times_ideal = wt_off_cell_ideal;
  waiting_time_struct(iter).name = bursting_sim_struct(s).name;
  iter = iter + 1;
end

% save
save([WritePath 'waiting_time_struct.mat'],'waiting_time_struct')