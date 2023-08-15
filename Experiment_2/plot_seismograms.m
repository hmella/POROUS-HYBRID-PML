close all

%% LOAD SEISMOGRAMS
% Simulations
simulations = {'extended';'paraxial';'hybrid';'hybrid_multiaxial'};
names = {'Extended';'Paraxial';'Hybrid PML';'Hybrid M-PML'};

% Number of seismograms
seismograms = [1,5];
Nb_seism = numel(seismograms);

% Cell object which contains the seismograms info
s = cell([numel(simulations) Nb_seism]);
m = cell([numel(simulations) Nb_seism]); % magnitudes

% Load seismograms of all simulations
for i=1:numel(simulations)
    for j=1:Nb_seism
        try
          s{i,j} = load(sprintf('energy and seismograms/s%d_%s.txt',seismograms(j),simulations{i}));
          m{i,j} = [s{i,j}(:,1),...
                    sqrt(s{i,j}(:,2).^2 + s{i,j}(:,3).^2),...
                    sqrt(s{i,j}(:,4).^2 + s{i,j}(:,5).^2),...
                    abs(s{i,j}(:,6))];
        catch
          s{i,j} = NaN([100, 6]);
          m{i,j} = NaN([100, 4]);
        end
    end
end
    
% Calculate differences
d = cell([numel(simulations) Nb_seism]); % differences
nb_elements = 1e+10*ones([1 numel(simulations)]);
for i=1:numel(simulations)
    for j=1:Nb_seism
        nb_elements(i) = numel(m{i,j}(:,1)); % nb of points in energy array
        if nb_elements(1) < nb_elements(i)
            range = 1:nb_elements(1);
        else
            range = 1:nb_elements(i);
        end
        % Normalized error
        d{i,j} = [s{i,j}(range,1),...
                sqrt((s{i,j}(range,2)-s{1,j}(range,2)).^2 + (s{i,j}(range,3)-s{1,j}(range,3)).^2)/max(m{1,j}(range,2)),...
                sqrt((s{i,j}(range,4)-s{1,j}(range,4)).^2 + (s{i,j}(range,5)-s{1,j}(range,5)).^2)/max(m{1,j}(range,3)),...
                abs(s{i,j}(range,6)-s{1,j}(range,6))/max(m{1,j}(range,4))];
     end
end


%% SEISMOGRAMS PLOTS
% Colors for plots
addpath('../../matlab_tools/linspecer/')
c = linspecer(5);

% Titles and linestyles for plots
variables = {'$|u|$ (m)';'$|w|$ (m)';'$p$ (Pa)'};
lines = {'-';'-.';'--';'--'};

% Plot options
linewidth = 1.3;
fontsize = 18;

% Letters for plots
letters = 'abcd';

% Plot seismograms
for seism=1:Nb_seism

    % Create figure
    figure(seism)
    tiledlayout(3,1,'Padding','compact','TileSpacing','compact');

    % Plot seismograms
    for var=1:numel(variables)

        % Plot paraxial
        nexttile(var)
        plot(m{1,seism}(:,1),m{1,seism}(:,var+1),lines{1},'LineWidth',linewidth,'Color',c(1,:)); hold on
        plot(m{2,seism}(:,1),m{2,seism}(:,var+1),lines{2},'LineWidth',linewidth,'Color',c(2,:)); hold on
        plot(m{3,seism}(:,1),m{3,seism}(:,var+1),lines{3},'LineWidth',linewidth,'Color',c(4,:)); hold off
        
        % Axis labels
        ylabel(variables{var},'Interpreter','latex','FontSize',fontsize)
        ax = gca;
        b_ylim = 1.1*max([m{1,seism}(:,var+1)',m{2,seism}(:,var+1)',m{3,seism}(:,var+1)']);
        a_ylim = -b_ylim/10;
        ax.XLim = [0 2.0];
        ax.YLim = [a_ylim b_ylim];
        
        % Legend and formating
        if var==3
            xlabel('Time (s)','Interpreter','latex','FontSize',fontsize)
            l = legend(names);
            l.Orientation = 'vertical';
            l.Interpreter = 'latex';
            l.Location = 'northeast';
            l.FontSize = fontsize;
        else
            ax.XAxis.TickLabels = [];
        end

        % Formating
        ax = gca;
        ax.XGrid = 'off';
        ax.YGrid = 'off';
        ax.Box = 'on';
        ax.FontWeight = 'bold';
        ax.FontSmoothing = 'on';
        ax.FontSize = fontsize;
        ax.LineWidth = 1.5;
        ax.TickLength = [0.025, 0.25];
        ax.XAxis.TickLength = [0.025, 0.25];
        ax.XAxis.TickValues = [0, 0.5, 1.0, 1.5, 2.0];
        ax.TickDir = 'in';
        ax.TickLabelInterpreter = 'latex';

    end

    % Formating and print
    set(ax,'Position',[0.072184483471957,0.101598796201987,0.880315516528043,0.220097200263579]);
    set(gcf,'Position',[577,282,911,663]);
    drawnow
    print('-depsc','-r300',sprintf('figures/Fig_8%s',letters(seism)))

end



%% ERROR PLOTS
% Titles and linestyles for plots
variables_error = {'$e_{u}$ (a.u.)';'$e_{w}$ (a.u.)';'$e_{p}$ (a.u.)'};
lines = {'-';'-';'-';'--'};

% Plot options
linewidth = 1.3;
fontsize = 18;

% Plot seismograms
for seism=1:Nb_seism

    % Create figure
    figure(Nb_seism+seism)
    tiledlayout(3,1,'Padding','compact','TileSpacing','compact');

    % Plot seismograms
    for var=1:numel(variables)

        % Plot errors
        nexttile(var)
        semilogy(d{2,seism}(:,1),d{2,seism}(:,var+1),lines{2},'LineWidth',linewidth,'Color',c(2,:)); hold on
        semilogy(d{3,seism}(:,1),d{3,seism}(:,var+1),lines{3},'LineWidth',linewidth,'Color',c(4,:)); hold on
        semilogy(d{4,seism}(:,1),d{4,seism}(:,var+1),lines{3},'LineWidth',linewidth,'Color',c(5,:)); hold off

        % Axis labels
        ylabel(variables_error{var},'Interpreter','latex','FontSize',fontsize)
        ax = gca;
        ax.XLim = [0 2.0];
        ax.YLim = [1e-8 5e2];
%         ax.YLim = [a_ylim b_ylim];

        % Legend and formating
        if var==3
            l = legend(names(2:4));
            l.Orientation = 'horizontal';
            l.Interpreter = 'latex';
            l.Location = 'northwest';
            l.FontSize = fontsize;
            xlabel('Time (s)','Interpreter','latex','FontSize',fontsize)
        else
            ax.XAxis.TickLabels = [];
        end

        % Formating
        ax = gca;
        ax.XGrid = 'off';
        ax.YGrid = 'off';
        ax.Box = 'on';
        ax.FontWeight = 'bold';
        ax.FontSmoothing = 'on';
        ax.FontSize = fontsize;
        ax.LineWidth = 1.5;
        ax.TickLength = [0.025, 0.25];
        ax.XAxis.TickLength = [0.025, 0.25];
        ax.XAxis.TickValues = [0, 0.5, 1.0, 1.5, 2.0];
        ax.YAxis.TickValues = [1e-10, 1e-5, 1];
        ax.TickDir = 'in';
        ax.TickLabelInterpreter = 'latex';

    end

    % Formating and print
    set(ax,'Position',[0.072184483471957,0.101598796201987,0.880315516528043,0.220097200263579]);
    set(gcf,'Position',[577,282,911,663]);
    drawnow
    print('-depsc','-r300',sprintf('figures/Fig_9%s',letters(seism)))
    

end