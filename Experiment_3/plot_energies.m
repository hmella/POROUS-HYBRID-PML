close all

%% Data import
% Simulations
simulations = {'extended';'paraxial';'mixed';'hybrid';'hybrid_multiaxial'};
names = {'Extended';'Paraxial';'Direct PML';'Hybrid PML';'Hybrid M-PML'};

% Arrays to store times and energies
times = cell([1, numel(simulations)]);
energies = cell(size(times));

% Load energies
for i=1:numel(simulations)
    try
        tmp = load(sprintf('energy and seismograms/energy_%s.txt',simulations{i}));
    catch
        tmp = NaN([10,2]);
    end
    times{i} = tmp(:,1);
    energies{i} = tmp(:,2);
end


%% PLOTS FOR EXTENDED, PARAXIAL, MIXED, AND HYBRID
% Make folder for plots
mkdir('figures/')

% Colors for plots
addpath('../../matlab_tools/linspecer/')
c = linspecer(5);

% Markers or lines for plots
linetypes = {'-','-.',':','--'};

% Line width
lwidth = 2.5;

% All methods
figure(1)
for i=1:numel(simulations(1:4))
    semilogy(times{i},energies{i},linetypes{i},'LineWidth',lwidth,'Color',c(i,:)); hold on
end
hold off
axis([-0.05 4.0 0 1000])
xlabel('Time (s)','Interpreter','latex','FontSize',14)
ylabel('Total energy (J)','Interpreter','latex','FontSize',14)

% Formating
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'off';
ax.Box = 'on';
ax.FontWeight = 'bold';
ax.FontSmoothing = 'on';
ax.FontSize = 14;
ax.LineWidth = 1.5;
ax.TickLength = [0.025, 0.25];
ax.XAxis.TickLength = [0.025, 0.25];
ax.TickDir = 'in';
ax.TickLabelInterpreter = 'latex';
l = legend(names);
l.Interpreter = 'latex';
l.Location = 'northeast';
l.FontSize = 14;
set(ax,'Position',[0.130000000000000,0.167764096665994,0.775000000000000,0.757235903334006]);
set(gcf,'Position',[675,632,733,325]);
drawnow
print('-depsc','-r300','figures/Fig_11a')


%% PLOTS FOR PML VARIATIONS
% Colors for plots
c = c(4:5,:);

% Markers or lines for plots
linetypes = {'-','-.',':','--'};

% Line width
lwidth = 2.5;

% PML plots
figure(2)
for i=1:numel(simulations(4:end))
    semilogy(times{i+3},energies{i+3},linetypes{i},'LineWidth',lwidth,'Color',c(i,:)); hold on
end
hold off
axis([-0.05 4.0 0 1e+3])
xlabel('Time (s)','Interpreter','latex','FontSize',14)
ylabel('Total energy (J)','Interpreter','latex','FontSize',14)

% Formating
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'off';
ax.Box = 'on';
ax.FontWeight = 'bold';
ax.FontSmoothing = 'on';
ax.FontSize = 14;
ax.LineWidth = 1.5;
ax.TickLength = [0.025, 0.25];
ax.XAxis.TickLength = [0.025, 0.25];
ax.TickDir = 'in';
ax.TickLabelInterpreter = 'latex';
l = legend(names(4:end));
l.Interpreter = 'latex';
l.Location = 'northeast';
l.FontSize = 14;
set(ax,'Position',[0.130000000000000,0.167764096665994,0.775000000000000,0.757235903334006]);
set(gcf,'Position',[675,632,733,325]);
drawnow
print('-depsc','-r300','figures/Fig_11b')

