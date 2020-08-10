% script to simulate burst time series with parameters from numerical and
% analytic results
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

% set basic plot parameters
t_max = 60;
ylimTrace = [-0.5 6.5];
n_bound_vec = 0:n_bcd_sites;

% sim name cell
sim_name_cell = {bursting_sim_struct.name};

% define colors
blue = [190 201 224]/255;
red = [246 141 100]/255;
green = [203 220 170]/255;
gray = [0.7020    0.7020    0.7020];
cmap1 = [green ; blue ;red];

% define resampling time res
resamp_freq = 1; % 1/seconds

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) independent binding (null model)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% specify appropriate index
ind_sim_index = find(contains(sim_name_cell,'independent'));

%  extract n bound vec
ind_plot_indices = 1:3;

% plot results of stochastic simulations
trace_index = 2;

state_fig = figure;
hold on
% generate resampled trace (moving average, essentially)
for i = 1%[2 1 3]%ind_plot_indices
  time_raw = double(bursting_sim_struct(ind_sim_index).sim_time_cell{i,trace_index});
  trace_raw = double(bursting_sim_struct(ind_sim_index).sim_emission_cell{i,trace_index});
  [trace_rs, time_rs] = resample(trace_raw,time_raw,resamp_freq,'pchip');

  plot(time_rs/60, trace_rs,'Color',green,'LineWidth',1.5);
end
ylim(ylimTrace)
xlim([0 t_max])
ylabel('transcription rate')
xlabel('time (minutes)')
box on
set(gca,'Fontsize',14,'YTick',n_bound_vec)
p = plot(0,0);
ax = gca;
ax.YColor = 'black';
ax.XColor = 'black';
StandardFigurePBoC(p,gca);
state_fig.InvertHardcopy = 'off';
saveas(state_fig,[FigurePath 'no-coop_trace.png'])
saveas(state_fig,[FigurePath 'no-coop_trace.pdf'])


hist_bins = linspace(0,1,35);

hist_fig = figure('Position',[100 100 256 512]);

hold on
% generate resampled trace (moving average, essentially)
for i = 1%ind_plot_indices
  time_raw = double(bursting_sim_struct(ind_sim_index).sim_time_cell{i,trace_index});
  trace_raw = double(bursting_sim_struct(ind_sim_index).sim_emission_cell{i,trace_index});
  
  dur_vec = diff([time_raw t_max*60]);
  state_vec = trace_raw'+1;
  states_u = unique(state_vec);
  stateSums = accumarray(state_vec,dur_vec');
  stateSums = stateSums(stateSums~=0);
  stateSums = stateSums/sum(stateSums);
  stateShares = zeros(size(n_bound_vec));
  stateShares(ismember(n_bound_vec,states_u-1)) = stateSums;

  barh(n_bound_vec,stateShares,1,'FaceColor',green,'FaceAlpha',0.8);
    
end

xlabel('probability')
box on

p = plot(0,0);
ax = gca;
ax.YColor = 'black';
ax.XColor = 'black';
xlim([0 0.55])
set(gca,'Fontsize',14,'xtick',0:.25:.5)
ylim([n_bound_vec(1)-0.5 n_bound_vec(end)+0.5])
StandardFigurePBoC(p,gca);
hist_fig.InvertHardcopy = 'off';
saveas(hist_fig,[FigurePath 'no-coop_hist.png'])
saveas(hist_fig,[FigurePath 'no-coop_hist.pdf'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (b) cooperative binding 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
% specify appropriate index
coop_sim_index = find(contains(sim_name_cell,'kon-mediated'));

%  extract n bound vec
coop_plot_index = 6;%size(bursting_sim_struct(2).SS,2)-2;

% get average on and off rates
n_states = 7;
SS = bursting_sim_struct(coop_sim_index).SS(:,coop_plot_index);
Q = bursting_sim_struct(coop_sim_index).Q(:,:,coop_plot_index)';
a = ones(n_states); m1 = tril(a,-1); m2 = tril(a,-2); m3 = triu(a,1); m4 = triu(a,2); m5 = ~~eye(n_states);
off_col = Q(m3&~m4);
on_col = Q(m1&~m2);
% on_rates = on_col./(6-n_bound_vec(1:end-1))'
% on_mean = SS(1:end-1)'*(on_rates.*on_col) / sum(SS(1:end-1).*on_col)
% off_freqs = 
% on_mean = SS(1:end-1)'*(on_col./(n_states-1-n_bound_vec(1:end-1)')) / sum(SS(1:end-1))

%
% plot results of stochastic simulations

trace_index = 23;

state_fig = figure;
cmap2 = brewermap(9,'Set2');
hold on
% generate resampled trace (moving average, essentially)
time_raw = double(bursting_sim_struct(coop_sim_index).sim_time_cell{coop_plot_index,trace_index});
trace_raw = double(bursting_sim_struct(coop_sim_index).sim_emission_cell{coop_plot_index,trace_index});

stairs(time_raw/60, trace_raw,'Color',blue,'LineWidth',1.5);

ylim(ylimTrace)
xlim([0 t_max])
ylabel('transcription rate')
xlabel('time (minutes)')
box on
set(gca,'Fontsize',14,'YTick',n_bound_vec)
p = plot(0,0);
% ax = gca;
% ax.YColor = 'black';
% ax.XColor = 'black';
StandardFigurePBoC(p,gca);
state_fig.InvertHardcopy = 'off';
saveas(state_fig,[FigurePath 'coop_trace.png'])
saveas(state_fig,[FigurePath 'coop_trace.pdf'])

% calculate fraction of time in eac state
dur_vec = diff([time_raw t_max*60]);
stateSums = accumarray(trace_raw'+1,dur_vec');
stateShares = stateSums/sum(stateSums);

hist_fig = figure('Position',[100 100 256 512]);
hold on
barh(n_bound_vec,stateShares,1,'FaceColor',blue);
  
xlabel('probability')
box on
p = plot(0,0);
% ax = gca;
% ax.YColor = 'black';
xlim([0 0.55])
% ax.XColor = 'black';
set(gca,'Fontsize',14,'xtick',0:.25:.5)
ylim([n_bound_vec(1)-0.5 n_bound_vec(end)+0.5])
StandardFigurePBoC(p,gca);

hist_fig.InvertHardcopy = 'off';
saveas(hist_fig,[FigurePath 'coop_hist.png'])
saveas(hist_fig,[FigurePath 'coop_hist.pdf'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) rate-limiting step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% specify appropriate index
rateLim_sim_index = find(contains(sim_name_cell,'2rate-limiting'));

%  extract n bound vec
rateLim_plot_index = coop_plot_index;%size(bursting_sim_struct(rateLim_sim_index).SS,2);

% plot results of stochastic simulations
trace_index = 5;

state_fig = figure;
hold on

% generate resampled trace (moving average, essentially)
time_raw = double(bursting_sim_struct(rateLim_sim_index).sim_time_cell{rateLim_plot_index,trace_index});
trace_raw = double(bursting_sim_struct(rateLim_sim_index).sim_emission_cell{rateLim_plot_index,trace_index});

stairs(time_raw/60, trace_raw,'Color',red,'LineWidth',1.5);

ylim(ylimTrace)
xlim([0 t_max])
ylabel('transcription rate')
xlabel('time (minutes)')
box on
set(gca,'Fontsize',14,'YTick',n_bound_vec)
p = plot(0,0);
ax = gca;
ax.YColor = 'black';
ax.XColor = 'black';
StandardFigurePBoC(p,gca);
state_fig.InvertHardcopy = 'off';
saveas(state_fig,[FigurePath 'rateLim_trace.png'])
saveas(state_fig,[FigurePath 'rateLim_trace.pdf'])

% calculate fraction of time in eac state
dur_vec = diff([time_raw t_max*60]);
stateSums = accumarray(trace_raw'+1,dur_vec');
stateShares = stateSums/sum(stateSums);

hist_fig = figure('Position',[100 100 256 512]);
hold on
barh(n_bound_vec,stateShares,1,'FaceColor',red);
  
xlabel('probability')
box on
p = plot(0,0);
ax = gca;
ax.YColor = 'black';
ax.XColor = 'black';
xlim([0 0.55])
ylim([n_bound_vec(1)-0.5 n_bound_vec(end)+0.5])
StandardFigurePBoC(p,gca);
set(gca,'Fontsize',14,'xtick',0:.25:.5)
hist_fig.InvertHardcopy = 'off';
saveas(hist_fig,[FigurePath 'rateLim_hist.png'])
saveas(hist_fig,[FigurePath 'rateLim_hist.pdf'])
