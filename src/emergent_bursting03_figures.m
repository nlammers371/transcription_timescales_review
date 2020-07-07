% script to simulate burst time series with parameters from numerical and
% analytic results
clear
close all
cd('C:\Users\nlamm\projects\transcription_timescales_review\src')
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
load([DataPath 'bursting_calc_struct.mat'])

% set basic plot parameters
t_max = 60;
ylTrace = [-0.01 1.01];
% define colors
blue = [ 0.5529    0.6275    0.7961];
red = [0.9882    0.5529    0.3843];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) independent binding (null model)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
% blue = .7*[27 117 188]/256 + .3*[1 1 1];

% specify appropriate index
ind_sim_index = 1;

%  extract n bound vec
n_bound_vec = bursting_calc_struct(ind_sim_index).n_bound_vec;
ind_plot_indices = 1:2:size(bursting_sim_struct(ind_sim_index).SS,2);

% plot probability landscape
ind_probs_fig = figure;
cmap1 = brewermap(length(ind_plot_indices)+2,'Greys');
hold on
iter = 1;
for e = ind_plot_indices
  prob_vec = bursting_sim_struct(ind_sim_index).SS(:,e);  
  area((n_bound_vec)/n_bcd_sites, prob_vec,'FaceColor',cmap1(iter+1,:),'FaceAlpha',0.7)
  iter = iter+1;
end
ylim([0 0.42])
ylabel('state likelihood')
xlabel('normalized transcription rate')
box on
set(gca,'Fontsize',14)

ax = gca;
ax.YColor = 'black';
ax.XColor = 'black';
ax.Color = [228,221,209]/255;
ind_probs_fig.InvertHardcopy = 'off';
% grid on
saveas(ind_probs_fig,[FigurePath 'no-coop_probs.png'])
saveas(ind_probs_fig,[FigurePath 'no-coop_probs.pdf'])


% plot results of stochastic simulations
% plot parameters
trace_index = 2;

state_fig = figure;
hold on
for p = 2
    stairs(bursting_sim_struct(ind_sim_index).sim_time_cell{ind_plot_indices(p),trace_index}/60,...
    (bursting_sim_struct(ind_sim_index).sim_state_cell{ind_plot_indices(p),trace_index}-1)/n_bcd_sites,...
    'Color',cmap1(p+1,:),'LineWidth',0.5);
%   patchline(bursting_sim_struct(ind_sim_index).sim_time_cell{p,plot_index}/60,...
%     (bursting_sim_struct(ind_sim_index).sim_state_cell{p,plot_index}-1)/n_bcd_sites,...
%     'EdgeColor',cmap1(p,:),'LineWidth',1,'edgealpha',.5);
end
ylim(ylTrace)
xlim([0 t_max])
ylabel('normalized transcription rate')
xlabel('time (minutes)')
box on
set(gca,'Fontsize',14,'YTick',0:0.2:1)
p = plot(0,0);
ax = gca;
ax.YColor = 'black';
ax.XColor = 'black';
StandardFigurePBoC(p,gca);
state_fig.InvertHardcopy = 'off';
saveas(state_fig,[FigurePath 'no-coop_trace.png'])
saveas(state_fig,[FigurePath 'no-coop_trace.pdf'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (b) cooperative binding 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all

% specify appropriate index
coop_sim_index = 2;
coop_plot_index = 4;
% plot probability landscape
coop_probs_fig = figure;
hold on
iter = 1;
% first plot ind for reference
prob_vec = bursting_sim_struct(ind_sim_index).SS(:,ind_plot_indices(2));  
area((n_bound_vec)/n_bcd_sites, prob_vec,'FaceColor',cmap1(2+1,:),'FaceAlpha',0.5)

% now with cooperativity
prob_vec = bursting_sim_struct(coop_sim_index).SS(:,coop_plot_index);  
area((n_bound_vec)/n_bcd_sites, prob_vec,'FaceColor',blue,'FaceAlpha',0.7,'EdgeAlpha',1)
iter = iter + 1;


ylabel('state likelihood')
xlabel('normalized transcription rate')
box on
set(gca,'Fontsize',14)

ax = gca;
ax.YColor = 'black';
ax.XColor = 'black';
ax.Color = [228,221,209]/255;
coop_probs_fig.InvertHardcopy = 'off';

saveas(coop_probs_fig,[FigurePath 'coop_probs.png'])
saveas(coop_probs_fig,[FigurePath 'coop_probs.pdf'])


% plot results of stochastic simulations
% plot parameters
trace_index = 2;

state_fig = figure;
hold on

stairs(bursting_sim_struct(coop_sim_index).sim_time_cell{coop_plot_index,trace_index}/60,...
  (bursting_sim_struct(coop_sim_index).sim_state_cell{coop_plot_index,trace_index}-1)/n_bcd_sites,...
  'Color',blue,'LineWidth',1);

ylim(ylTrace)
xlim([0 t_max])
ylabel('normalized transcription rate')
xlabel('time (minutes)')
box on
set(gca,'Fontsize',14,'YTick',0:0.2:1)
p = plot(0,0);
ax = gca;
ax.YColor = 'black';
ax.XColor = 'black';
StandardFigurePBoC(p,gca);
state_fig.InvertHardcopy = 'off';
saveas(state_fig,[FigurePath 'coop_trace.png'])
saveas(state_fig,[FigurePath 'coop_trace.pdf'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) rate-limiting step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% close all

% specify appropriate index
bottleneck_sim_index = 4;

% plot probability landscape
bottleneck_probs_fig = figure;
hold on

% first plot ind for reference
prob_vec = bursting_sim_struct(ind_sim_index).SS(:,ind_plot_indices(1));  
area((n_bound_vec)/n_bcd_sites, prob_vec,'FaceColor',cmap1(2,:),'FaceAlpha',0.7)
prob_vec = bursting_sim_struct(ind_sim_index).SS(:,ind_plot_indices(3));  
area((n_bound_vec)/n_bcd_sites, prob_vec,'FaceColor',cmap1(end-1,:),'FaceAlpha',0.7)

% now with rate-limiting step
prob_vec = bursting_sim_struct(bottleneck_sim_index).SS(:,coop_plot_index);  
prob_vec = prob_vec(1:length(n_bound_vec)) + prob_vec(length(n_bound_vec)+1:end);
area((n_bound_vec)/n_bcd_sites, prob_vec,'FaceColor',red,'FaceAlpha',1,'EdgeAlpha',1)

ylabel('state likelihood')
xlabel('normalized transcription rate')
box on
set(gca,'Fontsize',14)
ylim([0 0.42])
ax = gca;
ax.YColor = 'black';
ax.XColor = 'black';
ax.Color = [228,221,209]/255;
bottleneck_probs_fig.InvertHardcopy = 'off';
saveas(bottleneck_probs_fig,[FigurePath 'bottleneck_probs.png'])
saveas(bottleneck_probs_fig,[FigurePath 'bottleneck_probs.pdf'])


% plot results of stochastic simulations
% plot parameters
trace_index = 8;
state_fig = figure;

% calculate effective state vec
state_vec = mod(bursting_sim_struct(bottleneck_sim_index).sim_state_cell{coop_plot_index,trace_index}-1,n_bound_vec(end)+1);

stairs(bursting_sim_struct(bottleneck_sim_index).sim_time_cell{coop_plot_index,trace_index}/60,...
  state_vec/n_bcd_sites,  'Color',red,'LineWidth',0.5);

ylim(ylTrace)
xlim([0 t_max])
ylabel('normalized transcription rate')
xlabel('time (minutes)')
box on
set(gca,'Fontsize',14,'YTick',0:0.2:1)
ax = gca;
ax.YColor = 'black';
ax.XColor = 'black';
ax.Color = [228,221,209]/255;
state_fig.InvertHardcopy = 'off';
saveas(state_fig,[FigurePath 'coop_trace.png'])
saveas(state_fig,[FigurePath 'coop_trace.pdf'])
