% analyze_results.m
clear; clc; close all;

% Define the experiment parameters to analyze.
target_sample = 1;
states_list = 2:6;
seeds_list = 1:10;
base_output_dir = '..\results';

% Initialize arrays to store the best results for model selection.
max_likelihoods = zeros(length(states_list), 1);
best_seeds = zeros(length(states_list), 1);

% ---------------------------------------------------------
% Part 1: Model Structure Estimation (Likelihood Analysis).
% ---------------------------------------------------------

% Iterate through each state size to find the best model.
for i = 1:length(states_list)
    s = states_list(i);
    current_best_ll = -Inf;
    current_best_seed = -1;
    
    % Check all seeds to avoid local minima and find the maximum likelihood.
    for seed = seeds_list
        % Construct the path to the likelihood file.
        target_dir = sprintf('%s/sample%d/state_%d/seed_%d', ...
                             base_output_dir, target_sample, s, seed);
        filename = sprintf('%s/likelihood%d.txt', target_dir, target_sample);
        
        % Read the likelihood value if the file exists.
        if exist(filename, 'file')
            val = load(filename);
            if val > current_best_ll
                current_best_ll = val;
                current_best_seed = seed;
            end
        end
    end
    
    % Store the best results for this state size.
    max_likelihoods(i) = current_best_ll;
    best_seeds(i) = current_best_seed;
end

% Plot the relationship between the number of states and log-likelihood.
figure(1);
plot(states_list, max_likelihoods, '-o', 'LineWidth', 2);
grid on;
title(['Model Selection for Sample ', num2str(target_sample)]);
xlabel('Number of States');
ylabel('Max Log-Likelihood');
xticks(states_list);

% ---------------------------------------------------------
% Part 2: Convergence Analysis (Learning Curve).
% ---------------------------------------------------------

% Select the best model (e.g., State=4) to visualize its convergence.
% Here we automatically pick the state with the highest overall likelihood.
[~, best_idx] = max(max_likelihoods);
best_state = states_list(best_idx);
best_seed_overall = best_seeds(best_idx);

% Construct the path to the convergence history file of the best model.
conv_dir = sprintf('%s/sample%d/state_%d/seed_%d', ...
                   base_output_dir, target_sample, best_state, best_seed_overall);
conv_file = sprintf('%s/convergence%d.txt', conv_dir, target_sample);

% Load and plot the convergence data.
if exist(conv_file, 'file')
    data = load(conv_file); % Columns: [Iteration, Error]
    
    figure(2);
    plot(data(:,1), data(:,2), 'LineWidth', 1.5);
    grid on;
    title(sprintf('Convergence (State=%d, Seed=%d)', best_state, best_seed_overall));
    xlabel('Iteration');
    ylabel('Total Error');
else
    fprintf('Convergence file not found for the best model.\n');
end