% main_plot_figures.m
% Figure generation script for HMM experiments.
% This script generates all figures described in figure.txt.

clear; clc; close all;

%% Path configuration.
scriptDir = fileparts(mfilename('fullpath'));
rootDir   = fileparts(scriptDir);

resultsDir = fullfile(rootDir, 'results');
figDir     = fullfile(resultsDir, 'figures');
sampleName = 'sample0';

if ~exist(figDir, 'dir')
    mkdir(figDir);
end

masterPath = fullfile(resultsDir, 'master.tsv');

%% Load master table.
opts = detectImportOptions(masterPath, 'FileType', 'text');
opts = setvartype(opts, 'char');
T = readtable(masterPath, opts);

% Filter sample.
T = T(strcmp(T.SampleName, sampleName), :);

stateList = unique(T.States);

%% ================================
% Figure 1: Stability Analysis.
% ================================
figure;
hold on;

for k = stateList'
    idx = T.States == k;
    ll  = T.FinalLogLik(idx);
    ll  = sort(ll);
    plot(ll, 'LineWidth', 1.5);
end

xlabel('Trial Index (Sorted)');
ylabel('Final Log-Likelihood');
legend(arrayfun(@(x) sprintf('K=%d', x), stateList, 'UniformOutput', false), ...
       'Location', 'best');
grid on;
set(gca, 'Box', 'on');

saveas(gcf, fullfile(figDir, 'fig1_stability_sorted.png'));

%% ================================
% Figure 2: Learning Convergence.
% ================================
figure;
hold on;

for k = stateList'
    idxK = T.States == k;

    % Select best seed.
    [~, id] = max(T.FinalLogLik(idxK));
    row = T(idxK, :);
    histPath = row.HistoryPath{id};

    H = readtable(histPath, 'FileType', 'text');
    semilogy(H.Step, abs(H.DError), 'LineWidth', 1.5);
end

xlabel('Iteration');
ylabel('Log Error');
legend(arrayfun(@(x) sprintf('K=%d (Best)', x), stateList, 'UniformOutput', false), ...
       'Location', 'best');
grid on;
set(gca, 'Box', 'on');

saveas(gcf, fullfile(figDir, 'fig2_convergence.png'));

%% ================================
% Figure 3: Model Selection.
% ================================
figure;

maxLL  = zeros(numel(stateList), 1);
minAIC = zeros(numel(stateList), 1);
minBIC = zeros(numel(stateList), 1);

for i = 1:numel(stateList)
    k = stateList(i);
    idx = T.States == k;
    maxLL(i)  = max(T.FinalLogLik(idx));
    minAIC(i) = min(T.AIC(idx));
    minBIC(i) = min(T.BIC(idx));
end

yyaxis left;
plot(stateList, minAIC, '-o', 'LineWidth', 1.5); hold on;
plot(stateList, minBIC, '-s', 'LineWidth', 1.5);
ylabel('AIC / BIC');

yyaxis right;
plot(stateList, maxLL, '-^', 'LineWidth', 1.5);
ylabel('Max Log-Likelihood');

xlabel('Number of States');
legend({'AIC', 'BIC', 'Log-Likelihood'}, 'Location', 'best');
grid on;
set(gca, 'Box', 'on');

saveas(gcf, fullfile(figDir, 'fig3_model_selection.png'));

%% ================================
% Figure 4: Structure Analysis.
% ================================

% Determine optimal K by BIC.
[~, idxBest] = min(minBIC);
Kbest = stateList(idxBest);

idx = T.States == Kbest;
[~, id] = max(T.FinalLogLik(idx));
row = T(idx, :);

modelPath   = row.ModelPath{id};
viterbiPath = fullfile(resultsDir, sampleName, 'viterbi', ...
    sprintf('s%d_%d_out.txt', Kbest, row.Seed(id)));

%% Load HMM model.
fid = fopen(modelPath, 'r');
A = [];
B = [];
pi = [];

while ~feof(fid)
    line = strtrim(fgetl(fid));
    if startsWith(line, 'State')
        vals = split(line);
        pi(end+1) = str2double(vals{3});
        fgetl(fid); % Output
        outLine = strtrim(ans);
    end
    if startsWith(line, 'Output')
        B(end+1, :) = str2double(split(line(8:end)));
    end
    if startsWith(line, 'Transition')
        A(end+1, :) = str2double(split(line(12:end)));
    end
end
fclose(fid);

%% Load Viterbi path.
V = readmatrix(viterbiPath);

figure;

subplot(3,1,1);
imagesc(A);
colorbar;
title('Transition Matrix A');
xlabel('To State');
ylabel('From State');

subplot(3,1,2);
imagesc(B);
colorbar;
title('Emission Matrix B');
xlabel('Symbol');
ylabel('State');

subplot(3,1,3);
plot(V, 'LineWidth', 1.5);
xlabel('Time');
ylabel('State');
title('Viterbi Path');

saveas(gcf, fullfile(figDir, 'fig4_structure_analysis.png'));

fprintf('All figures generated successfully.\n');
