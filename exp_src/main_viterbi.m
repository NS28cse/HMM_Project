% main_viterbi.m
% Automation script for running Viterbi decoding only.
% This script assumes that best_model_s[States].txt already exists.

clear; clc; close all;

%% Experiment configuration.
samples   = {'sample0', 'sample1', 'sample2', 'sample3', 'sample4'};
stateList = 2:6;

% Seed is only used for naming consistency in Viterbi.
% It does NOT affect decoding if the model is fixed.
seed = 0;

%% Path configuration.
scriptDir = fileparts(mfilename('fullpath'));
rootDir   = fileparts(scriptDir);

srcDir     = fullfile(rootDir, 'src');
binDir     = fullfile(rootDir, 'bin');
dataDir    = fullfile(rootDir, 'data', 'raw', 'samples');
resultsDir = fullfile(rootDir, 'results');

exeViterbi = fullfile(binDir, 'Viterbi.exe');

%% Compile Viterbi source.
if ~exist(binDir, 'dir')
    mkdir(binDir);
end

cppFile = fullfile(srcDir, 'Viterbi.cpp');
exeFile = exeViterbi;

cmd = sprintf('g++ -O3 "%s" -o "%s"', cppFile, exeFile);
fprintf('Compiling %s.\n', cppFile);

[status, out] = system(cmd);
if status ~= 0
    error('Compilation failed:\n%s', out);
end

%% Viterbi execution loop.
for s = 1:numel(samples)
    sampleName = samples{s};
    fprintf('Processing sample: %s\n', sampleName);

    sampleDataDir = fullfile(dataDir, sampleName);
    sampleResDir  = fullfile(resultsDir, sampleName);

    modelDir   = fullfile(sampleResDir, 'models');
    viterbiDir = fullfile(sampleResDir, 'viterbi');

    if ~exist(viterbiDir, 'dir')
        mkdir(viterbiDir);
    end

    for K = stateList
        bestModelPath = fullfile(modelDir, sprintf('best_model_s%d.txt', K));

        if ~exist(bestModelPath, 'file')
            warning('Best model not found: %s', bestModelPath);
            continue;
        end

        fprintf('Running Viterbi: %s States=%d\n', sampleName, K);

        cmd = sprintf('"%s" "%s" "%s" "%s" %d %d', ...
            exeViterbi, sampleDataDir, viterbiDir, bestModelPath, K, seed);

        [status, out] = system(cmd);
        if status ~= 0
            warning('Viterbi failed: %s States=%d.', sampleName, K);
            continue;
        end

        files = dir(fullfile(viterbiDir, sprintf('s%d_*_out.txt', K)));
        fprintf('>Viterbi finished (%d sequences)\n', numel(files));
    end
end

fprintf('All Viterbi processes finished successfully.\n');
