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
DataPath = ['../out/waiting_time_distributions/' project '/'];
% DataPathSim = ['../out/emergent_bursting/' project '/'];

% load data
load([DataPath 'waiting_time_struct.mat'])

% set basic plot parameters
t_max = 60;
ylimTrace = [-0.5 6.5];
n_bound_vec = 0:n_bcd_sites;

% sim name cell
sim_name_cell = {waiting_time_struct.name};

% define colors
purple = [171 133 172]/255;
blue = [190 201 224]/255;
red = [246 141 100]/255;
green = [203 220 170]/255;
gray = [0.7020    0.7020    0.7020];
cmap1 = [green ; blue ;red];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1) Make figure illustrating passage time concept
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% specify appropriate index
rateLim_sim_index = find(contains(sim_name_cell,'2rate-limiting'));
%
%  extract n bound vec
rateLim_sim_sub_index = 2;
% plot results of stochastic simulations
trace_index = 6;

state_fig = figure;%('Position',[100 100 1024 512]);
hold on

% extract trace data
% time_raw = double(bursting_sim_struct(rateLim_sim_index).sim_time_cell{rateLim_plot_index,trace_index});
% trace_raw = double(bursting_sim_struct(rateLim_sim_index).sim_emission_cell{rateLim_plot_index,trace_index});
trace_raw = waiting_time_struct(rateLim_sim_index).trace_array(trace_index,:,rateLim_sim_sub_index)-1;

% extract corresponding viterbi fit
viterbi_time = waiting_time_struct(rateLim_sim_index).time_vector;
viterbi_fit = waiting_time_struct(rateLim_sim_index).viterbi_traces(trace_index,:,rateLim_sim_sub_index)*n_bcd_sites;

stairs(viterbi_time/60, trace_raw,'Color',purple,'LineWidth',1);
stairs(viterbi_time/60, viterbi_fit,'Color','k','LineWidth',1.5);

ylim(ylimTrace)
xlim([0 t_max])
ylabel('transcription rate')
xlabel('time (minutes)')
box on
% legend('raw trace','2 state fit')
set(gca,'Fontsize',14,'YTick',n_bound_vec)
p = plot(0,0);
ax = gca;
ax.YColor = 'black';
ax.XColor = 'black';
StandardFigurePBoC(p,gca);
state_fig.InvertHardcopy = 'off';
saveas(state_fig,[FigurePath 'rateLim_trace.png'])
saveas(state_fig,[FigurePath 'rateLim_trace.pdf'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (2) Make figures showing passage times for rate-limiting step mode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_step_vec = [1 2 5 15];
rateLim_sim_indices = find(ismember(sim_name_cell,{'1rate-limiting steps','2rate-limiting steps','5rate-limiting steps','15rate-limiting steps'}));
close all


% define bins for grouping waiting time measurements
wt_bins = linspace(0,4,50);
wt_centers = (wt_bins(1:end-1)+wt_bins(2:end))/2;

% define gamma function
gamma_fun = @(x) x(1) * x(2)^x(3) .* wt_centers.^(x(3)-1).*exp(-x(2)*wt_centers) ./ gamma(x(3));
exp_fun = @(x) x(1) * exp(-x(2) .* wt_centers);
options = optimoptions(@lsqnonlin,'MaxFunctionEvaluations',5e4);
% perform fits
rl_fit_struct = struct;
for i = 1:length(rateLim_sim_indices)
  wt_vec_raw = waiting_time_struct(rateLim_sim_indices(i)).off_waiting_times_ideal{1};
  wt_vec = wt_vec_raw / mean(wt_vec_raw);
  p_vec = histcounts(wt_vec,wt_bins);
  if i == 1
    fit_obj = @(x) exp_fun(x)-p_vec;
    rl_fit_struct(i).params = lsqnonlin(fit_obj,[1 1], [0 0],[Inf Inf],options);
    pd_vec = exp_fun(rl_fit_struct(i).params);
  else
    fit_obj = @(x) gamma_fun(x)-p_vec;
    rl_fit_struct(i).params = lsqnonlin(fit_obj,[1 1 1], [0 0 0],[Inf Inf Inf],options);
    pd_vec = gamma_fun(rl_fit_struct(i).params);
  end  
  rl_fit_struct(i).pd_vec = pd_vec/sum(pd_vec);
  rl_fit_struct(i).p_vec = p_vec / sum(p_vec);
  rl_fit_struct(i).var = var(wt_vec_raw);
  rl_fit_struct(i).mean = mean(wt_vec_raw);
end



offset = 2;
rl_hist_fig = figure;
cmap1_raw = brewermap([],'Pastel2');
cmap2_raw = brewermap([],'Set2');
start1 = brighten(cmap1_raw(2,:),-0.8);
end1 = 0.5*brighten(cmap1_raw(2,:),0.8)+[.5 .5 .5];
start2 = brighten(cmap2_raw(2,:),-.8);
end2 = .5*brighten(cmap2_raw(2,:),.8)+[.5 .5 .5];


cmap1 = interp1([1,length(rateLim_sim_indices)],[start1 ; end1],linspace(1,length(rateLim_sim_indices),length(rateLim_sim_indices)+offset));
cmap2 = interp1([1,length(rateLim_sim_indices)],[start2 ; end2],linspace(1,length(rateLim_sim_indices),length(rateLim_sim_indices)+offset));

hold on
iter = length(rateLim_sim_indices);
pl = [];
lgd_str = {};
for i = fliplr(1:length(rateLim_sim_indices))
  p_vec = rl_fit_struct(i).p_vec;  
  bar(wt_centers,p_vec,1,'FaceColor',cmap1(i+offset/2,:),'FaceAlpha',.7,'EdgeAlpha',0.3);   
  pd_vec = rl_fit_struct(i).pd_vec;
  pl = [pl plot(wt_centers,pd_vec,'Color',cmap2(i+offset/2,:),'LineWidth',2)]; 
  if i == 1
    lgd_str = [lgd_str{:} {[num2str(n_step_vec(i)) ' rate-limiting step']}];
  else
    lgd_str = [lgd_str{:} {[num2str(n_step_vec(i)) ' rate-limiting steps']}];
  end
end
% for i = 1:length(rateLim_sim_indices)
%   pd_vec = rl_fit_struct(i).pd_vec;
%   plot(wt_centers,pd_vec,'Color',cmap(i,:),'LineWidth',2);  
% end

xlim([0 wt_bins(end)])
ylim([0 .13])
set(gca,'Fontsize',14)
ylabel('share')
xlabel('normalized passage time from 0->6')
p = plot(0,0);
legend(pl,lgd_str{:});
ax = gca;
ax.YColor = 'black';
ax.XColor = 'black';
StandardFigurePBoC(p,gca);
rl_hist_fig.InvertHardcopy = 'off';
saveas(rl_hist_fig,[FigurePath 'rateLim_wt_hist.png'])
saveas(rl_hist_fig,[FigurePath 'rateLim_wt_hist.pdf'])



coop_sim_index = find(ismember(sim_name_cell,{'koff-mediated cooperativity'}));
sim_vec = 1:length(waiting_time_struct(i).off_waiting_times_ideal);
plot_index = length(waiting_time_struct(i).off_waiting_times_ideal);

blue = [190 201 224]/255;
blueDark = brighten(blue,-0.5);

% perform fits
coop_fit_struct = struct;
for i = sim_vec
  wt_vec_raw = waiting_time_struct(coop_sim_index).off_waiting_times_ideal{i};
  wt_vec = wt_vec_raw / mean(wt_vec_raw );    
  p_vec = histcounts(wt_vec,wt_bins);
  
  fit_obj = @(x) exp_fun(x)-p_vec;
  coop_fit_struct(i).params = lsqnonlin(fit_obj,[1 1], [0 0],[Inf Inf],options);
  pd_vec = exp_fun(coop_fit_struct(i).params);

  coop_fit_struct(i).pd_vec = pd_vec/sum(pd_vec);
  coop_fit_struct(i).p_vec = p_vec / sum(p_vec);
  coop_fit_struct(i).var = var(wt_vec_raw);
  coop_fit_struct(i).mean = mean(wt_vec_raw);
end


coop_hist_fig = figure;
offset = length(rateLim_sim_indices);
hold on

p_vec = coop_fit_struct(plot_index).p_vec;  
bar(wt_centers,p_vec,1,'FaceColor',blue,'FaceAlpha',.7,'EdgeAlpha',0.3);   
pd_vec = coop_fit_struct(plot_index).pd_vec;
plot(wt_centers,pd_vec,'Color',blueDark,'LineWidth',2);  

xlim([0 wt_bins(end)])
set(gca,'Fontsize',14)
ylabel('share')
xlabel('normalized passage time from 0->6')
p = plot(0,0);
ax = gca;
ax.YColor = 'black';
ax.XColor = 'black';
StandardFigurePBoC(p,gca);
coop_hist_fig.InvertHardcopy = 'off';
saveas(coop_hist_fig,[FigurePath 'coop_wt_hist.png'])
saveas(coop_hist_fig,[FigurePath 'coop_wt_hist.pdf'])


%% make fano factor figure



