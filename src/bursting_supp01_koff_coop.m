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
FigurePath = ['../fig/bursting_supp/' project '/'];
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

% define resampling time res. Slower sampling time is neeeded to make plots
% intelligible
resamp_res = 0.5; % in seconds

% define time grid for resampling
time_rs = 0:resamp_res:3600;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% koff-mediated cooperative binding 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% specify appropriate index
coop_sim_index = find(contains(sim_name_cell,'koff-mediated'));

coop_plot_index = 1; 
sim_param_indices_coop = 176:4:201;
% plot results of stochastic simulations
trace_index = 16; 

state_fig = figure;
cmap2 = brewermap(9,'Set2');
hold on

% generate resampled trace (moving average, essentially)
time_raw = repelem(double(bursting_sim_struct(coop_sim_index).sim_time_cell{coop_plot_index,trace_index}),1);
trace_raw = repelem(double(bursting_sim_struct(coop_sim_index).sim_emission_cell{coop_plot_index,trace_index}),1);

trace_rs = interp1(time_raw,trace_raw,time_rs,'previous');

stairs(time_rs/60, trace_rs,'Color',[blue 0.0],'LineWidth',1.5);

ylim(ylimTrace)
xlim([0 t_max])
ylabel('transcription rate')
xlabel('time (minutes)')
box on
set(gca,'Fontsize',14,'YTick',n_bound_vec)
p = plot(0,0);
StandardFigurePBoC(p,gca);
state_fig.InvertHardcopy = 'off';

saveas(state_fig,[FigurePath 'coop_trace.png'])
saveas(state_fig,[FigurePath 'coop_trace.pdf'])