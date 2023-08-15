clear; clc; close all

% Experiments folder
experiments = {'Experiment_1/';'Experiment_2/';'Experiment_3/'};

% Run matlab scripts
for i=1:numel(experiments)
    cd(experiments{i})
    plot_energies
    plot_seismograms
    cd ..
end