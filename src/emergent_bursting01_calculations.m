% script to solve for effeCooptive off/on rates in linear chain binding model
clear
close all
cd('C:\Users\nlamm\projects\transcription_timescales_review\src')
addpath('utilities')

% Since we are assuming that our system is a linear markov chain, it follows 
% that the chain must be in thermodynamic equilbrium. % Thus we will take 
% the approach of first calculating how the relative occupancies of each 
% state in the chain change as a function of binding/unbinding cooparativity 
% terms. Kinetics can then be speCoopified by invoking measured timescales for
% Bcd unbinding (Mir et al 2018)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define basic hyperparameters and symbols
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_bs = 6; % number of Bcd binding sites
n_states = n_bs + 1; % total number of states
n_med = 3; % midpoint 
n_bound_vec = 0:n_bs; % vector encoding # bound in each state
n_bound_vec_plot = linspace(0,n_bs);
off_rate_basal = 1/3; % this sets overall system timescales (from Mir et al, 2018)
n_calc_points = 101; % number of ref points
simIndices = [1 round(prctile(1:n_calc_points,[10 50 90])) n_calc_points];
% set save paths
project = ['n' num2str(n_bs) ]; % projeCoopt identifier

DataPath = ['../out/emergent_bursting/' project '/'];
mkdir(DataPath)

%%% Calculate transition rate matrices and state probabilities for:
%%% (a) independent binding (null model)
%%% (b) cooperative binding
%%% (c) rate-limiting step

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% account for state multplicities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Since our states correspond to a certain # bound and NOT the
% bound/unbound state of each site, there are multiple ways to realize each
% state. For instance, there are 6! / (4! x 2!) = 15 ways to get 2 Bcd moleCoopules
% bound to the enhancer. These multiplicities must be accounted for in our
% effeCooptive state energies

mult_vec = NaN(size(n_bound_vec));
for w = 1:length(mult_vec)
  mult_vec(w) = factorial(n_bs) ./ (factorial(n_bs-n_bound_vec(w)) .*factorial(n_bound_vec(w)));
end
mu_vec = -log(mult_vec);

% initialize energy variables. "eBoundBcd" indicates the energy associated with a 
% single bound bcd moleCoopule. "eCoop" indicates the amount of
% cooperativity/synergy that exists between Bcd molecules 
syms eBoundBcd eCoop

% define function to calculate state energies without cooperativity
ind_state_energy_vec = matlabFunction(n_bound_vec*eBoundBcd + mu_vec);

% plot energy landscape for a variety of eBoundBcd values

eBoundBcd_vals = linspace(-2,2,n_calc_points); % different binding energies to explore

% perform all these calculations assuming eBcd = 0
ebIndex = ceil(n_calc_points/2);

% define function to calculate state energies with cooperativity. Note that
% for convenience, we define the cooperative contribution such that the
% energy potential has the form of a quadratic peaked at n_med. This is an
% arbitrary choice, but represents perhaps the simplest energy landscape
% that is capable of generating bimodal behaviors necessary for bursting
coop_state_energy_vec = matlabFunction(eCoop*(n_bound_vec-n_med).^2);% - ...
%   realmin*(eBoundBcd*n_bound_vec - mu_vec));
eCoop_vals = -linspace(0,1,n_calc_points); % different cooperativity energies to explore

%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate state probabilities and on/off rates 
%%%%%%%%%%%%%%%%%%%%%%%%%%

% (a) for independent case, basal off rate and binding affinity alone
%     determine the transition kinetics between states

% (b) There are two different "flavors" of cooperativity: off rate-mediated
%     cooeprativity wherein bound factors stabilize each-other and therefore
%     decrease the effective unbinding rate, and on rate-mediated
%     cooperativity, wherein having factors bound somehow increases the
%     propensity for additional factors to bind.

% (c) Here the rate of the slow step sets overall kinetics. binding
%     energies within open and closed states dictate kinetics therein


% for each energy landscape, we can calculate the implied on rate for a 
% single Bcd molecule (absent cooperative effects)
on_rate_basal_vec = NaN(1,n_calc_points); 

% initialize arrays to store transition rate (generator) matrices for each
% condition
q_ind_array = zeros(n_bs+1,n_bs+1,n_calc_points); % no cooperativity
q_coop_off_array = zeros(n_bs+1,n_bs+1,n_calc_points); % off rate-mediated
q_coop_on_array = zeros(n_bs+1,n_bs+1,n_calc_points); % on rate-mediated

% individual state probabilities
state_prob_ind_array = NaN(length(n_bound_vec),n_calc_points);
state_prob_coop_array = NaN(length(n_bound_vec),n_calc_points);

% define matrices for indexing
a = ones(n_bs+1);
m1 = tril(a,-1);
m2 = tril(a,-2);
m3 = triu(a,1);
m4 = triu(a,2);
m5 = ~~eye(n_bs+1);

% calculate metrics for (a) and (b) first
for e = 1:length(eBoundBcd_vals)
  % no cooperativity
  state_probs_ind = exp(-ind_state_energy_vec(eBoundBcd_vals(e)));
  state_probs_ind = state_probs_ind/sum(state_probs_ind);
  state_prob_ind_array(:,e) = state_probs_ind;
  
  % cooperativity
  state_probs_coop = exp(-(ind_state_energy_vec(eBoundBcd_vals(ebIndex))+coop_state_energy_vec(eCoop_vals(e))));
  state_probs_coop = state_probs_coop/sum(state_probs_coop);
  state_prob_coop_array(:,e) = state_probs_coop;
  
  % calculate basal on rate using detailed balance
  on_rate_basal_vec(e) = off_rate_basal * state_probs_ind(2) / state_probs_ind(1) / n_bs;
  
  % state-speCoopific rates (independent)
  rate_ind_temp(:,1) = off_rate_basal * n_bound_vec;
  rate_ind_temp(:,2) = on_rate_basal_vec(e) * (n_bs-n_bound_vec);
  
  % generator matrix for ind
  ind_slice = zeros(n_bs+1);
  ind_slice(m1&~m2) = rate_ind_temp(1:end-1,2);
  ind_slice(m3&~m4) = rate_ind_temp(2:end,1);
  ind_slice(m5) = -sum(ind_slice);
  q_ind_array(:,:,e) = ind_slice';
  
  % calculate state-specific on and off rates for on and off rate mediated
  % cooperativities
  
  % on rate-mediated
  rate_coop_on_temp(:,1) = off_rate_basal * n_bound_vec;
  rate_coop_on_temp(1:end-1,2) = rate_coop_on_temp(2:end,1)' .* state_probs_coop(2:end) ./ state_probs_coop(1:end-1);
  
  % generator matrix for ind
  c_on_slice = zeros(n_bs+1);
  c_on_slice(m1&~m2) = rate_coop_on_temp(1:end-1,2);
  c_on_slice(m3&~m4) = rate_coop_on_temp(2:end,1);
  c_on_slice(m5) = -sum(c_on_slice);
  q_coop_on_array(:,:,e) = c_on_slice';
  
  % off rate-mediated
  rate_coop_off_temp(:,2) = on_rate_basal_vec(e) * (n_bs-n_bound_vec);
  rate_coop_off_temp(2:end,1) = rate_coop_off_temp(1:end-1,2)' .* state_probs_coop(1:end-1) ./ state_probs_coop(2:end);  
  
  % generator matrix for ind
  c_off_slice = zeros(n_bs+1);
  c_off_slice(m1&~m2) = rate_coop_off_temp(1:end-1,2);
  c_off_slice(m3&~m4) = rate_coop_off_temp(2:end,1);
  c_off_slice(m5) = -sum(c_off_slice);
  q_coop_off_array(:,:,e) = c_off_slice';
end


% solve for effective on and off rates. This is not strictly necessary but
% will be used as a consistency check that stochastic simulations are
% behaving as expected

% initialize arrays
eff_ton_ind_vec = NaN(1,n_calc_points);
eff_toff_ind_vec = NaN(1,n_calc_points);

eff_ton_on_coop_vec = NaN(1,n_calc_points);
eff_toff_on_coop_vec = NaN(1,n_calc_points);

eff_ton_off_coop_vec = NaN(1,n_calc_points);
eff_toff_off_coop_vec = NaN(1,n_calc_points);

% specify states from which to calculate on/off rates. This choice is
% somewhat arbitrary, but we know kon must be from a low state to a high
% state and koff must by high->low
calc_vec = [1 length(n_bound_vec)];

for e = 1:length(eBoundBcd_vals)
  % calculate effeCooptive off rates (5->3)
  
  % independent
  [eff_ton_ind_vec(e), eff_toff_ind_vec(e)] = pt_solve(q_ind_array(:,:,e),calc_vec(1),calc_vec(2));
  
  % koff-mediated cooperativity
  [eff_ton_on_coop_vec(e), eff_toff_on_coop_vec(e)] = pt_solve(q_coop_on_array(:,:,e),calc_vec(1),calc_vec(2));
  
  % kon-mediated cooperativity
  [eff_ton_off_coop_vec(e), eff_toff_off_coop_vec(e)] = pt_solve(q_coop_on_array(:,:,e),calc_vec(1),calc_vec(2)); 
  
end


% Now calculate metrics for (c)
eBcd_state_vec = [simIndices(1) simIndices(end)]; % define the two states to lie at the extrema

% set range of switching kinetics to explore
slow_kinetics_array = [1./eff_ton_on_coop_vec' 1./eff_ton_off_coop_vec']; 

% initialize array to store rate matrics
q_coop_rate_lim_array = zeros(2*(n_bs+1),2*(n_bs+1),n_calc_points); % on rate-mediated

% array to store individual state probabilities
state_prob_rate_lim_array = NaN(2*length(n_bound_vec),n_calc_points);

% calculate indices that correspond to exchange points 
sub1 = 1:length(n_bound_vec);
sub2 = length(n_bound_vec)+1:2*length(n_bound_vec);
linear_indices1 = sub2ind([2*(n_bs+1),2*(n_bs+1)],sub2,sub1);
linear_indices2 = sub2ind([2*(n_bs+1),2*(n_bs+1)],sub1,sub2);

for e = 1:n_calc_points
  % extract component rate arrays
  q1 = q_ind_array(:,:,eBcd_state_vec(1));
  q2 = q_ind_array(:,:,eBcd_state_vec(2));
  
  % initialize combined rate array
  q_combined = q_coop_rate_lim_array(:,:,e);
  q_combined(sub1,sub1) = q1';
  q_combined(sub2,sub2) = q2';
  % add exchange terms
  q_combined(linear_indices1) = slow_kinetics_array(e,1); % on rate
  q_combined(linear_indices2) = slow_kinetics_array(e,2); % off rate
  q_combined(eye(size(q_combined))==1) = 0;
  q_combined(eye(size(q_combined))==1) = -sum(q_combined);
  q_coop_rate_lim_array(:,:,e) = q_combined';
  % calculate state occupancies
  [V,D] = eig(q_combined);
  [~,mi] = max(diag(D));
  ss_vec = V(:,mi)/sum(V(:,mi));
  state_prob_rate_lim_array(:,e) = ss_vec;%(sub1) + ss_vec(sub2);
end

%%% save data structures to use for simulations
bursting_calc_struct = struct;

%(a) independent binding/unbinding 
bursting_calc_struct(1).name = 'independent';
bursting_calc_struct(1).state_prob_array = state_prob_ind_array;
bursting_calc_struct(1).rate_array = q_ind_array;
bursting_calc_struct(1).eff_on_states = calc_vec;
bursting_calc_struct(1).eff_on_rates = [1./eff_ton_ind_vec' 1./eff_toff_ind_vec'];

%(b.1) kon-mediated binding/unbinding 
bursting_calc_struct(2).name = 'kon-mediated cooperativity';
bursting_calc_struct(2).state_prob_array = state_prob_coop_array;
bursting_calc_struct(2).rate_array = q_coop_on_array;
bursting_calc_struct(2).eff_on_states = calc_vec;
bursting_calc_struct(2).eff_rates = [1./eff_ton_on_coop_vec' 1./eff_toff_on_coop_vec'];

%(b.2) kon-mediated binding/unbinding 
bursting_calc_struct(3).name = 'koff-mediated cooperativity';
bursting_calc_struct(3).state_prob_array = state_prob_coop_array;
bursting_calc_struct(3).rate_array = q_coop_off_array;
bursting_calc_struct(3).eff_on_states = calc_vec;
bursting_calc_struct(3).eff_rates = [1./eff_ton_off_coop_vec' 1./eff_toff_off_coop_vec'];

%(c) rate-limiting step 
bursting_calc_struct(4).name = 'rate-limiting step';
bursting_calc_struct(4).state_prob_array = state_prob_rate_lim_array;
bursting_calc_struct(4).rate_array = q_coop_rate_lim_array;
bursting_calc_struct(4).eff_on_states = calc_vec;
bursting_calc_struct(4).eff_rates = [slow_kinetics_array' slow_kinetics_array']; % Note that this is only approximate (assumes strong separation of timescales)

for i = 1:4
  bursting_calc_struct(i).off_rate_basal = off_rate_basal;
  bursting_calc_struct(i).on_rate_basal_vec = on_rate_basal_vec;
  bursting_calc_struct(i).eCoop_vals = eCoop_vals;
  bursting_calc_struct(i).eBoundBcd_vals = eBoundBcd_vals;
  bursting_calc_struct(i).n_bound_vec = n_bound_vec;
  bursting_calc_struct(i).n_med = n_med;
  bursting_calc_struct(i).coop_state_energy_vec = coop_state_energy_vec;
  bursting_calc_struct(i).ind_state_energy_vec = ind_state_energy_vec;
  bursting_calc_struct(i).ebIndex = ebIndex;
  bursting_calc_struct(i).simIndices = simIndices;
end

save([DataPath, 'bursting_calc_struct.mat'],'bursting_calc_struct');
