% script to solve for effective off/on rates in linear chain binding model
clear
close all
cd('C:\Users\nlamm\projects\transcription_timescales_review\src')
addpath('utilities')

% Since we are assuming that our system is a linear markov chain, it follows that the chain must be in thermodynamic equilbrium.
% Thus we will take the approach of first calculating how the relative occupancies of each state in the chain change as a function
% of binding/unbinding cooparativity terms. Kinetics can then be specified by invoking known(ish) timescales for Bcd binding and unbinding 

% define basic hperparameters and symbols
n_bs = 6;
n_states = n_bs + 1;
n_med = 3;
n_bound_vec = 0:n_bs;
% account for state multplicities
w_vec = NaN(size(n_bound_vec));
for w = 1:length(w_vec)
  w_vec(w) = factorial(n_bs) ./ (factorial(n_bs-n_bound_vec(w)) .*factorial(n_bound_vec(w)));
end
mu_vec = -log(w_vec);
syms eb ec
%%
% specify simulation type
project = ['n' num2str(n_bs) ];
% set figure path
FigPath = ['../fig/emergent_bursting/' project '/'];
mkdir(FigPath)
DataPath = ['../out/emergent_bursting/' project '/'];
mkdir(DataPath)

% define state energies without cooperativity
ind_state_energy_vec = matlabFunction(n_bound_vec*eb + mu_vec);

% plot energy landscape for a variety of eb values
n_plot = 10;
eb_vals = linspace(-1,1,n_plot);

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% No cooperativity
%%%%%%%%%%%%%%%%%%%%%%%%%%

ind_fig = figure;
cmap = brewermap(n_plot,'Spectral');
hold on
for e = 1:length(eb_vals)
  plot(n_bound_vec, ind_state_energy_vec(eb_vals(e)) ,'Color',cmap(e,:))
  scatter(n_bound_vec, ind_state_energy_vec(eb_vals(e)),'MarkerEdgeAlpha',0,'MarkerFaceColor',cmap(e,:))
end
xlim([-.5 n_bs+.5])
ylabel('energy (au)')
xlabel('n bound')
box on
set(gca,'Fontsize',14)

ax = gca;
ax.YColor = 'black';
ax.XColor = 'black';
ax.Color = [228,221,209]/255;

ind_fig.InvertHardcopy = 'off';
saveas(ind_fig,[FigPath 'no-coop_energies.png'])
saveas(ind_fig,[FigPath 'no-coop_energies.pdf'])

ind_prob_fig = figure;
cmap = brewermap(n_plot,'Spectral');
hold on
for e = 1:length(eb_vals)
  prob_vec = exp(-ind_state_energy_vec(eb_vals(e)));
  prob_vec = prob_vec / sum(prob_vec);
  
  plot(n_bound_vec,  prob_vec,'Color',cmap(e,:))
  scatter(n_bound_vec, prob_vec ,'MarkerEdgeAlpha',0,'MarkerFaceColor',cmap(e,:))
end
xlim([-.5 n_bs+.5])
ylabel('state probabilities')
xlabel('n bound')
box on
set(gca,'Fontsize',14)

ax = gca;
ax.YColor = 'black';
ax.XColor = 'black';
ax.Color = [228,221,209]/255;
set(gca,'yScale','log')
ind_prob_fig.InvertHardcopy = 'off';
grid on
saveas(ind_prob_fig,[FigPath 'no-coop_probss.png'])
saveas(ind_prob_fig,[FigPath 'no-coop_probss.pdf'])

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% With cooperativity
%%%%%%%%%%%%%%%%%%%%%%%%%%

syms ec
close all
ec_vals = -linspace(0,1,n_plot);
coop_state_energy_vec = matlabFunction(ec*(n_bound_vec-n_med).^2 - eb*n_bound_vec - mu_vec);

coop_energy_fig = figure;
cmap = brewermap(n_plot,'Spectral');
hold on
for e = 1:length(eb_vals)
  state_energies = ind_state_energy_vec(eb_vals(e))+coop_state_energy_vec(eb_vals(e),ec_vals(e));
  plot(n_bound_vec,state_energies ,'Color',cmap(e,:))
  scatter(n_bound_vec,state_energies ,'MarkerEdgeAlpha',0,'MarkerFaceColor',cmap(e,:))
end
ylim([-10 1])
xlim([-.5 n_bs+.5])
ylabel('state energies (au)')
xlabel('n bound')
box on
set(gca,'Fontsize',14)
p = plot(0,0);
% set(gca,'yScale','log')
ax = gca;
ax.YColor = 'black';
ax.XColor = 'black';
grid on
ax.Color = [228,221,209]/255;
coop_energy_fig.InvertHardcopy = 'off';

saveas(coop_energy_fig,[FigPath 'coop_energies.png'])
saveas(coop_energy_fig,[FigPath 'coop_energies.pdf'])


coop_prob_fig = figure;
cmap = brewermap(n_plot,'Spectral');
hold on
for e = 1:length(eb_vals)
  state_probs_ind = exp(-(ind_state_energy_vec(eb_vals(e))+coop_state_energy_vec(eb_vals(e),ec_vals(e))));
  plot(n_bound_vec,state_probs_ind/sum(state_probs_ind),'Color',cmap(e,:))
  scatter(n_bound_vec,state_probs_ind/sum(state_probs_ind) ,'MarkerEdgeAlpha',0,'MarkerFaceColor',cmap(e,:))
end
% ylim([-10 1])
xlim([-.5 n_bs+.5])
ylabel('state probabilities')
xlabel('n bound')
box on
set(gca,'Fontsize',14)
p = plot(0,0);
set(gca,'yScale','log')
ax = gca;
ax.YColor = 'black';
ax.XColor = 'black';
grid on
ax.Color = [228,221,209]/255;
coop_prob_fig.InvertHardcopy = 'off';

saveas(coop_prob_fig,[FigPath 'coop_probs.png'])
saveas(coop_prob_fig,[FigPath 'coop_probs.pdf'])

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% estimate on and off rates assuming (1) on-rate and (2) off-rate meidated cooperativity
%%%%%%%%%%%%%%%%%%%%%%%%%%

off_rate_basal = 1/3; % this sets overall system timescales

% transition rate vectors
on_rate_basal_vec = NaN(1,n_plot);
rate_ind_array = zeros(length(n_bound_vec),2,n_plot);
q_ind_array = zeros(n_bs+1,n_bs+1);
rate_coop_off_array = zeros(length(n_bound_vec),2,n_plot);
q_coop_off_array = zeros(n_bs+1,n_bs+1);
rate_coop_on_array = zeros(length(n_bound_vec),2,n_plot);
q_coop_on_array = zeros(n_bs+1,n_bs+1);

% individual state probabilities
state_prob_ind_array = NaN(length(n_bound_vec),n_plot);
state_prob_coop_array = NaN(length(n_bound_vec),n_plot);

% mixing rate vector
mix_rate_ind_vec = NaN(1,n_plot);
mix_rate_coop_vec = NaN(1,n_plot);

% define matrices for indexing
a = ones(n_bs+1);
m1 = tril(a,-1);
m2 = tril(a,-2);
m3 = triu(a,1);
m4 = triu(a,2);
m5 = ~~eye(n_bs+1);

for e = 1:length(eb_vals)
  % no cooperativity
  state_probs_ind = exp(-ind_state_energy_vec(eb_vals(e)));
  state_probs_ind = state_probs_ind/sum(state_probs_ind);
  state_prob_ind_array(:,e) = state_probs_ind;
  
  % cooperativity
  state_probs_coop = exp(-(ind_state_energy_vec(eb_vals(e))+coop_state_energy_vec(eb_vals(e),ec_vals(e))));
  state_probs_coop = state_probs_coop/sum(state_probs_coop);
  state_prob_coop_array(:,e) = state_probs_coop;
  
  % calculate basal on rate using detailed balance
  on_rate_basal_vec(e) = off_rate_basal * state_probs_ind(2) / state_probs_ind(1) / n_bs;
  
  % state-specific rates (independent)
  rate_ind_array(:,1,e) = off_rate_basal * n_bound_vec;
  rate_ind_array(:,2,e) = on_rate_basal_vec(e) * (n_bs-n_bound_vec);
  
  % generator matrix for ind
  ind_slice = zeros(n_bs+1);
  ind_slice(m1&~m2) = rate_ind_array(1:end-1,2,e);
  ind_slice(m3&~m4) = rate_ind_array(2:end,1,e);
  ind_slice(m5) = -sum(ind_slice);
  q_ind_array(:,:,e) = ind_slice;
  
  % calculate state-specific on and off rates for on and off rate mediated
  % cooperativities
  
  % on rate-mediated
  rate_coop_off_array(:,1,e) = off_rate_basal * n_bound_vec;
  rate_coop_off_array(1:end-1,2,e) = rate_coop_off_array(2:end,1,e)' .* state_probs_coop(2:end) ./ state_probs_coop(1:end-1);
  
  % generator matrix for ind
  c_off_slice = zeros(n_bs+1);
  c_off_slice(m1&~m2) = rate_coop_off_array(1:end-1,2,e);
  c_off_slice(m3&~m4) = rate_coop_off_array(2:end,1,e);
  c_off_slice(m5) = -sum(c_off_slice);
  q_coop_off_array(:,:,e) = c_off_slice;
  
  % off rate-mediated
  rate_coop_on_array(:,2,e) = on_rate_basal_vec(e) * (n_bs-n_bound_vec);
  rate_coop_on_array(2:end,1,e) = rate_coop_on_array(1:end-1,2,e)' .* state_probs_coop(1:end-1) ./ state_probs_coop(2:end);  
  
  % generator matrix for ind
  c_on_slice = zeros(n_bs+1);
  c_on_slice(m1&~m2) = rate_coop_on_array(1:end-1,2,e);
  c_on_slice(m3&~m4) = rate_coop_on_array(2:end,1,e);
  c_on_slice(m5) = -sum(c_on_slice);
  q_coop_on_array(:,:,e) = c_on_slice;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% now solve for first passage times
%%%%%%%%%%%%%%%%%%%%%%%%%%
% define function to calculate effective off rates


% initialize arrays
eff_kon_ind_vec = NaN(1,n_plot);
eff_koff_ind_vec = NaN(1,n_plot);
eff_kon_on_coop_vec = NaN(1,n_plot);
eff_koff_on_coop_vec = NaN(1,n_plot);
eff_kon_off_coop_vec = NaN(1,n_plot);
eff_koff_off_coop_vec = NaN(1,n_plot);

% define symbolic vectors
syms ton [n_bs+1 1]
ton(n_med+2) = 0;

syms toff [n_bs+1 1]
toff(n_med) = 0;
tic
for e = 1:length(eb_vals)
  % calculate effective off rates (5->3)
  
  % independent
  q_ind = q_ind_array(:,:,e)';
  eff_kon_ind_vec(e) = 1/pt_solve(q_ind,n_med+2,n_med);
  eff_koff_ind_vec(e) = 1/pt_solve(q_ind,n_med,n_med+2);
  
  % koff coop
  q_koff = q_coop_off_array(:,:,e)';
  eff_kon_off_coop_vec(e) = 1/pt_solve(q_koff,n_med+2,n_med);
  eff_koff_off_coop_vec(e) = 1/pt_solve(q_koff,n_med,n_med+2);
  
  % kon coop
  q_kon = q_coop_on_array(:,:,e)';
  eff_kon_on_coop_vec(e) = 1/pt_solve(q_kon,n_med+2,n_med);
  eff_koff_on_coop_vec(e) = 1/pt_solve(q_kon,n_med,n_med+2);
  
end

toc














function out = pt_solve(q_in,to,from)
  dim = size(q_in,1);
  syms t [dim 1]
  t(to) = 0;
  q = -diag(q_in);
  p = q_in ./ q;
  p(eye(dim)==1) = 0;
  p(to,:) = 0;
  p(:,to) = 0;
  q(to) = 0;
  sys_ind = p*t + q - t;
  pt_sol = solve(sys_ind,t([1:to-1 to+1:dim]));
  eval(['out = double(pt_sol.t' num2str(from) ');'])
end