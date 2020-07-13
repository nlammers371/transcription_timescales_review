% script to make figures examine waiting time distributions for 
% rate-limitiung step and and emergent cooperativity simulations
clear
close all
addpath('utilities')

% load numeric results
n_bcd_sites = 6;
project = ['n' num2str(n_bcd_sites)];
addpath('utilities')

% set paths
FigurePath = ['../fig/waiting_time_distributions/' project '/'];
mkdir(FigurePath)
DataPathWT = ['../out/waiting_time_distributions/' project '/'];
DataPathSim = ['../out/emergent_bursting/' project '/'];

% load data
load([DataPathWT 'waiting_time_struct.mat'])
load([DataPathSim 'bursting_sim_struct.mat'])

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (1) Make figure illustrating passage time concept
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% specify appropriate index
rateLim_sim_index = find(contains(sim_name_cell,'2rate-limiting'));

%  extract n bound vec
rateLim_plot_index = size(bursting_sim_struct(2).SS,2);

% plot results of stochastic simulations
trace_index = 5;

state_fig = figure;
hold on

% extract trace data
% time_raw = double(bursting_sim_struct(rateLim_sim_index).sim_time_cell{rateLim_plot_index,trace_index});
% trace_raw = double(bursting_sim_struct(rateLim_sim_index).sim_emission_cell{rateLim_plot_index,trace_index});
trace_raw = waiting_time_struct.trace_array(trace_index,:,rateLim_sim_index)-1;

% extract corresponding viterbi fit
viterbi_time = waiting_time_struct.time_vector;
viterbi_fit = waiting_time_struct.viterbi_traces(trace_index,:,rateLim_sim_index)*n_bcd_sites;

stairs(viterbi_time/60, trace_raw,'Color',red,'LineWidth',1);
stairs(viterbi_time/60, viterbi_fit,'Color','k','LineWidth',1);

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
% saveas(state_fig,[FigurePath 'rateLim_trace.png'])
% saveas(state_fig,[FigurePath 'rateLim_trace.pdf'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (2) Make figures showing passage times for different realizations of the two models
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
