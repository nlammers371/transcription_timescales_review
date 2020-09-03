% script to simulate full binding/unbinding model to check average
% unbinding rate
clear
close all
addpath('utilities')

% load numeric results
n_bcd_sites = 6;
project = ['n' num2str(n_bcd_sites)];
addpath('utilities')

% set paths
FigurePath = ['../fig/bursting_supp/' project '/'];
mkdir(FigurePath)
DataPath = ['../out/emergent_bursting/' project '/'];

% load data
load([DataPath 'bursting_sim_struct.mat'])
load([DataPath 'bursting_chain_calc_struct.mat'])


% sim name cell
sim_name_cell = {bursting_sim_struct.name};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% koff-mediated cooperative binding 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% specify appropriate index
coop_sim_index = find(contains(sim_name_cell,'koff-mediated'));

coop_plot_index = 1; 
sim_param_indices_coop = 176:4:201;

% extract transition rate matrix
Q = bursting_sim_struct(coop_sim_index).Q(:,:,coop_plot_index)';
SS = bursting_sim_struct(coop_sim_index).SS(:,coop_plot_index)';

% extract microscopic rates
n_vec = 0:n_bcd_sites;
n_states = length(n_vec);
a = ones(n_states); m1 = tril(a,-1); m2 = tril(a,-2); m3 = triu(a,1); m4 = triu(a,2); m5 = ~~eye(n_states);

% effective rates
k_minus_vec = [0 Q(m3&~m4)'];
k_plus_vec = [Q(m1&~m2)' 0];

% extract microscopic rates
k_unbind_vec = k_minus_vec ./n_vec;
k_unbind_vec(1) = 0;
k_bind_vec = k_plus_vec ./(n_bcd_sites-n_vec);
k_bind_vec(end) = 0;

%% build model
n_sim = 10; % number of simulations
T = 1e4; % total time to simulate in seconds
% cell arrays to store results
event_time_cell = cell(1,n_sim);
bound_state_cell = cell(1,n_sim);

for n = 1:n_sim
  % initialize
  event_time_vec = [0];
  
  nBound = randsample(n_vec,1,true,SS);
  s_bound = randsample(1:n_bcd_sites,nBound,false);
  bound_state_array = zeros(1,n_bcd_sites);
  bound_state_array(s_bound) = 1;
  
  simTime = 0; 
  
  while simTime < T
    
    currState = bound_state_array(end,:);
    nBound = sum(currState);
    
    % draw expected jump time
    tau = 1 / (k_minus_vec(nBound+1) + k_plus_vec(nBound+1));
    dt = exprnd(tau);
    
    simTime = simTime + dt;
    
    % determine whether we bind or unbind
    eventType = randsample([-1 1],1,true,[k_minus_vec(nBound+1) k_plus_vec(nBound+1)]);
    
    if eventType == 1      
      options = find(currState==0);      
    else
      options = find(currState==1);
    end
    if isempty(options)
      error('asda')
    end
    % select binding site to update    
    bs = randsample(repelem(options,2),1,false);% repelem prevents undesirable behavior when only 1 option
    currState(bs) = currState(bs) + eventType;
    if any(currState>1) || any(currState<0)
      error('afa')
    end
    bound_state_array(end+1,:) = currState;
    event_time_vec(end+1) = simTime;    
  end
  
  bound_state_cell{n} = bound_state_array;
  event_time_cell{n} = event_time_vec;
  
end

%% extract single-molecules unbinding times
unbinding_vec = [];
binding_vec = [];
for n = 1:n_sim
  bound_state_array = bound_state_cell{n};
  event_time_vec = event_time_cell{n};
  
  for i = 1:n_bcd_sites
    binding_state_vec = bound_state_array(:,i);
    
    change_vec = [0 diff(binding_state_vec')];
    
    binding_times_raw = event_time_vec(change_vec==1);
    unbinding_times_raw = event_time_vec(change_vec==-1);
     
    if ~isempty(binding_times_raw) && ~isempty(unbinding_times_raw)            
      binding_times1 = binding_times_raw(binding_times_raw<unbinding_times_raw(end));
      unbinding_times1 = unbinding_times_raw(unbinding_times_raw>binding_times_raw(1));
      
      unbinding_vec = [unbinding_vec unbinding_times1-binding_times1];
      
      binding_times2 = binding_times_raw(binding_times_raw>unbinding_times_raw(1));
      unbinding_times2 = unbinding_times_raw(unbinding_times_raw<binding_times_raw(end));
      
      binding_vec = [binding_vec binding_times2-unbinding_times2];
    end
  end   
    
end

%% fit exponential 
ub_bins = linspace(0,max(unbinding_vec),5e4);
ub_centers = ub_bins(1:end-1)+diff(ub_bins)/2;

p_vec = histcounts(unbinding_vec,ub_bins);
% p_vec = p_vec/sum(p_vec);

exp_fun = @(x) x(1) * exp(-ub_centers./x(2));
options = optimoptions(@lsqnonlin,'MaxFunctionEvaluations',5e4);
fit_obj = @(x) exp_fun(x)-p_vec;

params = lsqnonlin(fit_obj,[1 1], [0 0],[Inf Inf],options);

%% Make plots
% downsample to make this plotable
ds_vec = linspace(ub_centers(1),3,1e2);
p_vec_plot = interp1(ub_centers,p_vec,ds_vec); 
% re-adjust sum
p_vec_plot = p_vec_plot*sum(p_vec(ub_centers<=3))/sum(p_vec_plot);
close all

% make exponential fit figure first
bar_fig = figure;
cmap = brewermap(9,'Set2');

hold on
bar(ds_vec,p_vec_plot/sum(p_vec),1,'FaceColor',cmap(3,:),'FaceAlpha',1,'EdgeAlpha',1);   

pd_vec = exp_fun(params);
plot(ub_centers,pd_vec/sum(pd_vec),'Color','k','LineWidth',2);

xlim([0 2])
legend('activator dwell times','exponential fit (\mu = 0.23s)')
xlabel('activator dwell times (s)')
ylabel('probability')
% box on
set(gca,'Fontsize',14)
StandardFigurePBoC([],gca);
bar_fig.InvertHardcopy = 'off';
saveas(bar_fig,[FigurePath 'activator_dwell_exp.png'])
saveas(bar_fig,[FigurePath 'activator_dwell_exp.pdf'],'pdf')



%% now try to capture non-exponential behavior
log_bins = logspace(-3,log10(max(unbinding_vec)),3e2);
log_centers = log_bins(1:end-1)+diff(log_bins)/2;
log_counts = histcounts(unbinding_vec,log_bins);

pwr_fig = figure;
cmap = brewermap(9,'Set2');
hold on
% scatter(log_centers,log_counts)
b = bar(log_counts/sum(log_counts),1,'FaceColor',cmap(3,:),'FaceAlpha',1,'EdgeAlpha',1); 
set(gca,'Yscale','log')

% find where to place tick marks
xTickLabelsStr = {'10^{-3}'  '10^{-1}' '10^1' '10^3'};
xTickLabels = [10^-3  10^-1 10^1  10^3];
xTickPoints = [];
for x = 1:length(xTickLabels)
  [~,xTickPoints(x)] = min(abs(log_centers-xTickLabels(x)));
end

set(gca,'Xtick',xTickPoints,'xTickLabels',xTickLabelsStr)

xlabel('activator dwell times (sec)')
ylabel('probability')

box off
set(gca,'Fontsize',14)
StandardFigurePBoC([],gca);
pwr_fig.InvertHardcopy = 'off';

saveas(pwr_fig,[FigurePath 'activator_dwell_pwr.png'])
saveas(pwr_fig,[FigurePath 'activator_dwell_pwr.pdf'])