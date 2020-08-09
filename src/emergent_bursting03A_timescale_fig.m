% script to make burst model calculations appendix figure
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
load([DataPath 'bursting_chain_calc_struct.mat'])
load([DataPath 'bursting_step_calc_struct.mat'])

% define colors
pboc = [228,221,209]/255;
purple = brighten([171 133 172]/255,.5);
blue = [190 201 224]/255;
red = [246 141 100]/255;
green = [203 220 170]/255;
cmap0 = [green ; blue ;red];

% key calc quantities
n_bound_vec = 0:6;

% specify index to show
coop_index = 2;

% cooperativity energy
ec_vec = -bursting_chain_calc_struct(coop_index).activatorEnergies * 2/(n_bcd_sites-1);
omega_vec = exp(-ec_vec);
ec_filter = ec_vec <=0; % select for only synergistic interactions

% define quantities to plot
konEff = bursting_chain_calc_struct(coop_index).eff_rates(:,1);
koffEff = bursting_chain_calc_struct(coop_index).eff_rates(:,2);
kappa_vec =  (konEff+koffEff)./(konEff.*koffEff) / 60;

end_fraction = sum(bursting_chain_calc_struct(coop_index).SS([1 end],:));


et_ss_fig = figure;
hold on
cmap2 = brewermap(9,'Set2');

yyaxis right
% plot(omega_vec(ec_filter),konEff(ec_filter),'Color',cmap2(2,:),'LineWidth',1.5)
plot(omega_vec(ec_filter),kappa_vec(ec_filter),'Color',cmap2(3,:),'LineWidth',2)
% set(gca,'yScale','log','YColor','k')
set(gca,'YColor','k')
ylabel('bursting timescale (minutes)')
ylim([0.2 20])

yyaxis left
plot(omega_vec(ec_filter),end_fraction(ec_filter),'Color','k','LineWidth',2)
ylabel('fraction of time in end states')
% ylim(ylimTrace)
xlim([1 7])
xlabel('cooperativity strength (\omega)')

box on
set(gca,'Fontsize',14)

ax = gca;
ax.YColor = 'black';
ax.XColor = 'black';
ax.Color = pboc;
% StandardFigurePBoC(p,gca);
et_ss_fig.InvertHardcopy = 'off';
saveas(et_ss_fig,[FigurePath 'et_vs_ss_vs_omega.png'])
saveas(et_ss_fig,[FigurePath 'et_vs_ss_vs_omega.pdf'])


