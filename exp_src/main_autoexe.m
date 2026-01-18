% run_experiments.m
clear; clc;

% Define the configuration for the experiment.
% Set the target sample number, range of states to test, and range of random seeds.
target_sample = 1;
states_list = 2:6;
seeds_list = 1:10;

% Define the paths for the executables and the base output directory.
% Note that paths are relative to the 'analysis' folder.
exe_bw = '..\bin\BaumWelch.exe';
exe_vt = '..\bin\Viterbi.exe';
base_output_dir = '..\results';

% Loop through each state size to test structural differences.
for s = states_list
    
    % Loop through each random seed to account for initialization dependency.
    for seed = seeds_list
        
        % Construct the specific output directory path for this condition.
        % Structure: ../results/sampleX/state_Y/seed_Z/
        current_dir = sprintf('%s/sample%d/state_%d/seed_%d', ...
                              base_output_dir, target_sample, s, seed);
        
        % Create the directory if it does not exist.
        if ~exist(current_dir, 'dir')
            mkdir(current_dir);
        end
        
        % Construct the command string for Baum-Welch training.
        % Arguments: [sample_num] [num_states] [output_folder] [seed]
        cmd_bw = sprintf('%s %d %d "%s" %d', ...
                         exe_bw, target_sample, s, current_dir, seed);
        
        % Construct the command string for Viterbi inference.
        % Arguments: [sample_num] [output_folder]
        cmd_vt = sprintf('%s %d "%s"', ...
                         exe_vt, target_sample, current_dir);
        
        % Execute the training process and check for errors.
        fprintf('Processing: State=%d, Seed=%d ... ', s, seed);
        [status_bw, ~] = system(cmd_bw);
        
        % Execute Viterbi only if training was successful.
        if status_bw == 0
            system(cmd_vt);
            fprintf('Done.\n');
        else
            fprintf('Failed during training.\n');
        end
        
    end
end