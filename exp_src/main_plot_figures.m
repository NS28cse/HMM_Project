% main_plot_figures.m
% Plot all figures required for the report (figure2.txt).
% Target sample: sample0.
% Figures are saved to results/figures/.
%
% MATLAB version: 2025a.

clear; clc; close all;

%% Configuration.
sampleName = 'sample4';
stateList  = 2:6;
resultDir  = fullfile('..', 'results');
figureDir  = fullfile(resultDir, sampleName, 'figures');

if ~exist(figureDir, 'dir')
    mkdir(figureDir);
end

fprintf('Generating figures for %s...\n', sampleName);

%% Load master table.
masterPath = fullfile(resultDir, 'master.tsv');
if ~exist(masterPath, 'file')
    error('Master file not found: %s', masterPath);
end

T = readtable(masterPath, 'FileType', 'text', 'Delimiter', '\t');

% Filter target sample.
T = T(strcmp(T.SampleName, sampleName), :);

if isempty(T)
    error('No data found for sample: %s', sampleName);
end

%% ------------------------------------------------------------
% Figure 1: Sorted final log-likelihood (stability analysis).
% Shows existence of local minima and justifies best seed selection.
% -------------------------------------------------------------
fprintf('  Generating Figure 1: Stability Analysis...\n');

figure('Position', [100, 100, 800, 600], 'Color', 'w');
hold on;

set(gcf, 'Color', 'w');          % Figure background.
ax = gca;
set(ax, ...
    'Color', 'w', ...            % Axes background.
    'XColor', 'k', ...           % X-axis (ticks, labels).
    'YColor', 'k', ...           % Y-axis.
    'LineWidth', 1.2, ...
    'FontSize', 11);

% Define line styles and markers for monochrome plotting.
lineStyles = {'-', '--', '-.', ':', '-'};
markers = {'o', 's', '+', 'x', 'h'};

for i = 1:length(stateList)
    K = stateList(i);
    idx = T.States == K;
    loglik = T.FinalLogLik(idx);
    loglik = sort(loglik, 'ascend');
    
    % Plot with distinct line style and marker for each K.
    plot(1:length(loglik), loglik, 'k', ...
         'LineStyle', lineStyles{1}, ...
         'LineWidth', 1.5, ...
         'Marker', markers{i}, ...
         'MarkerSize', 5, ...
         'MarkerFaceColor', 'k', ...
         'DisplayName', sprintf('K=%d', K));
end

hold off;

xlabel('Trial Index Sorted by Log-Likelihood', 'FontSize', 12);
ylabel('Final Log-Likelihood', 'FontSize', 12);
lgd = legend('Location', 'best', 'FontSize', 10);
set(lgd, 'Color', 'w', 'TextColor', 'k', 'EdgeColor', 'k');
grid on;
box on;

exportgraphics(gcf, fullfile(figureDir, [sampleName, '_fig1_stability.png']));
savefig(gcf, fullfile(figureDir, [sampleName, '_fig1_stability.fig']));
exportgraphics(gcf, fullfile(figureDir, [sampleName, '_fig1_stability.pdf']));

%% ------------------------------------------------------------
% Figure 2: Learning convergence curves (best seed only).
% Demonstrates that Baum-Welch algorithm converges correctly.
% -------------------------------------------------------------
fprintf('  Generating Figure 2: Learning Convergence...\n');

figure('Position', [100, 100, 800, 600], 'Color', 'w');
hold on;

set(gcf, 'Color', 'w');          % Figure background.
ax = gca;
set(ax, ...
    'Color', 'w', ...            % Axes background.
    'XColor', 'k', ...           % X-axis (ticks, labels).
    'YColor', 'k', ...           % Y-axis.
    'LineWidth', 1.2, ...
    'FontSize', 11);

for i = 1:length(stateList)
    K = stateList(i);
    
    % Select best seed for each state size (maximum log-likelihood).
    idxK = T.States == K;
    subT = T(idxK, :);
    [~, imax] = max(subT.FinalLogLik);
    bestSeed = subT.Seed(imax);
    
    % Load convergence history.
    historyPath = fullfile(resultDir, sampleName, 'history', ...
        sprintf('history_s%d_%d.tsv', K, bestSeed));
    
    if ~exist(historyPath, 'file')
        warning('History file not found: %s', historyPath);
        continue;
    end
    
    H = readtable(historyPath, 'FileType', 'text', 'Delimiter', '\t');
    
    % Plot error on semi-log scale.
    errors = abs(H.DError);
    errors(errors < 1e-300) = 1e-300;   % Numerical floor.
    semilogy(H.Step, errors, 'k', ...
             'LineStyle', lineStyles{1}, ...
             'LineWidth', 1.5, ...
             'Marker', markers{i}, ...
             'MarkerSize', 5, ...
             'MarkerFaceColor', 'k', ...
             'DisplayName', sprintf('K=%d', K));

    set(gca, 'XScale', 'log');
end

hold off;

xlabel('Iteration', 'FontSize', 12);
ylabel('Delta Log Likelihood Error', 'FontSize', 12);
% legend('Location', 'northeast', 'FontSize', 10);
lgd = legend('Location', 'best', 'FontSize', 10);
set(lgd, 'Color', 'w', 'TextColor', 'k', 'EdgeColor', 'k');

grid on;
box on;

exportgraphics(gcf, fullfile(figureDir, [sampleName, '_fig2_convergence.png']));
savefig(gcf, fullfile(figureDir, [sampleName, '_fig2_convergence.fig']));
exportgraphics(gcf, fullfile(figureDir, [sampleName, '_fig2_convergence.pdf']));

%% ------------------------------------------------------------
% Figure 3: Model selection (LogLik vs AIC/BIC).
% Demonstrates optimal state number selection based on information criteria.
% -------------------------------------------------------------
fprintf('  Generating Figure 3: Model Selection...\n');

figure('Position', [100, 100, 800, 600], 'Color', 'w');

set(gcf, 'Color', 'w');          % Figure background.
ax = gca;
set(ax, ...
    'Color', 'w', ...            % Axes background.
    'XColor', 'k', ...           % X-axis (ticks, labels).
    'YColor', 'k', ...           % Y-axis.
    'LineWidth', 1.2, ...
    'FontSize', 11);

% Calculate best metrics for each state.
maxLogLik = arrayfun(@(k) max(T.FinalLogLik(T.States==k)), stateList);
minAIC = arrayfun(@(k) min(T.AIC(T.States==k)), stateList);
minBIC = arrayfun(@(k) min(T.BIC(T.States==k)), stateList);

% Left axis: AIC and BIC.
yyaxis left;
plot(stateList, minAIC, 'k-o', ...
     'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'k', ...
     'DisplayName', 'AIC');
hold on;
plot(stateList, minBIC, 'k-s', ...
     'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'k', ...
     'DisplayName', 'BIC');
ylabel('AIC / BIC', 'FontSize', 12);
ax = gca;
ax.YAxis(1).Color = 'k';

% Right axis: Max log-likelihood.
yyaxis right;
plot(stateList, maxLogLik, 'k--h', ...
     'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'w', ...
     'DisplayName', 'Max Log-Likelihood');
ylabel('Max Log-Likelihood', 'FontSize', 12);
ax.YAxis(2).Color = 'k';

hold off;

xlabel('Number of States', 'FontSize', 12);
xticks(stateList);
lgd = legend('AIC', 'BIC', 'Max Log-Likelihood', 'Location', 'best', 'FontSize', 10);
set(lgd, 'Color', 'w', 'TextColor', 'k', 'EdgeColor', 'k');

grid on;
box on;

exportgraphics(gcf, fullfile(figureDir, [sampleName, '_fig3_model_selection.png']));
savefig(gcf, fullfile(figureDir, [sampleName, '_fig3_model_selection.fig']));
exportgraphics(gcf, fullfile(figureDir, [sampleName, '_fig3_model_selection.pdf']));


%% ------------------------------------------------------------
% Figures 4 & 5: Structure analysis of the optimal model.
% Shows transition and emission matrices of the best model.
% Optimal model is selected based on minimum BIC (avoids overfitting).
% -------------------------------------------------------------
fprintf('  Generating Figures 4 & 5: Structure Analysis...\n');

% Determine globally optimal model based on minimum BIC.
[~, ibest] = min(T.AIC);
bestK     = T.States(ibest);
bestSeed  = T.Seed(ibest);
bestAIC   = T.AIC(ibest);

fprintf('    Optimal model: K=%d, Seed=%d, BIC=%.2f\n', bestK, bestSeed, bestAIC);

% Locate model file.
modelPath = fullfile(resultDir, sampleName, 'models', ...
    sprintf('markov_output_s%d_%d.txt', bestK, bestSeed));

if ~exist(modelPath, 'file')
    error('Model file not found: %s', modelPath);
end

% Load HMM model.
[A, B, pi] = load_hmm_model(modelPath);
K = size(A, 1);
M = size(B, 2);

% --- Figure 4: Transition matrix A.
fprintf('    Saving Figure 4: Transition Matrix A...\n');

figure('Position', [100, 100, 700, 600], 'Color', 'w');

imagesc(A);
colormap(flipud(gray));
colorbar;
caxis([0, 1]);

set(gcf, 'Color', 'w');          % Figure background.
ax = gca;
set(ax, ...
    'Color', 'w', ...            % Axes background.
    'XColor', 'k', ...           % X-axis (ticks, labels).
    'YColor', 'k', ...           % Y-axis.
    'LineWidth', 1.2, ...
    'FontSize', 11);

cb = colorbar;

set(cb, ...
    'Color', 'k', ...        % Label color.
    'XColor', 'k', ...       % Tick color (horizontal).
    'YColor', 'k', ...       % Tick color (vertical).
    'FontSize', 11, ...
    'FontWeight', 'bold', ...
    'LineWidth', 1.2);


% Add text annotations.
for i = 1:K
    for j = 1:K
        val = A(i,j);
        if val > 0.5
            txtColor = 'w';   % Dark cell -> white text.
        else
            txtColor = 'k';   % Light cell -> black text.
        end
        text(j, i, sprintf('%.2f', val), ...
     'HorizontalAlignment', 'center', ...
     'Color', txtColor, ...
     'FontSize', 10, 'FontWeight', 'bold');

    end
end

% Configure axes.
stateLabels = cellstr(char('A' + (0:K-1)'));
xticks(1:K);
yticks(1:K);
xticklabels(stateLabels);
yticklabels(stateLabels);
xlabel('To State', 'FontSize', 12);
ylabel('From State', 'FontSize', 12);
axis square;

exportgraphics(gcf, fullfile(figureDir, [sampleName, '_fig4_matrix_A.png']));
savefig(gcf, fullfile(figureDir, [sampleName, '_fig4_matrix_A.fig']));
exportgraphics(gcf, fullfile(figureDir, [sampleName, '_fig4_matrix_A.pdf']));


% --- Figure 5: Emission matrix B.
fprintf('    Saving Figure 5: Emission Matrix B...\n');

figure('Position', [100, 100, 900, 600], 'Color', 'w');

imagesc(B);
colormap(flipud(gray));
colorbar;
caxis([0, 1]);

set(gcf, 'Color', 'w');          % Figure background.
ax = gca;
set(ax, ...
    'Color', 'w', ...            % Axes background.
    'XColor', 'k', ...           % X-axis (ticks, labels).
    'YColor', 'k', ...           % Y-axis.
    'LineWidth', 1.2, ...
    'FontSize', 11);

cb = colorbar;

set(cb, ...
    'Color', 'k', ...        % Label color.
    'XColor', 'k', ...       % Tick color (horizontal).
    'YColor', 'k', ...       % Tick color (vertical).
    'FontSize', 11, ...
    'FontWeight', 'bold', ...
    'LineWidth', 1.2);

% Add text annotations.
for i = 1:K
    for j = 1:M
        val = B(i,j);
        if val > 0.5
            txtColor = 'w';   % Dark cell -> white text.
        else
            txtColor = 'k';   % Light cell -> black text.
        end

        text(j, i, sprintf('%.2f', val), ...
        'HorizontalAlignment', 'center', ...
        'Color', txtColor, ...
        'FontSize', 10, 'FontWeight', 'bold');
    end
end

% Configure axes.
xticks(1:M);
yticks(1:K);
xticklabels(0:M-1);  % Symbols are 0-indexed.
yticklabels(stateLabels);
xlabel('Output Symbol', 'FontSize', 12);
ylabel('State', 'FontSize', 12);

exportgraphics(gcf, fullfile(figureDir, [sampleName, '_fig5_matrix_B.png']));
savefig(gcf, fullfile(figureDir, [sampleName, '_fig5_matrix_B.fig']));
exportgraphics(gcf, fullfile(figureDir, [sampleName, '_fig5_matrix_B.pdf']));

fprintf('All figures generated successfully.\n');
fprintf('Output directory: %s\n', figureDir);

%% ============================================================
% Local functions.
% ============================================================

function [A, B, pi] = load_hmm_model(filename)
% Load HMM parameters from a model text file.
% 
% Input:
%   filename - Path to the HMM model file.
% 
% Output:
%   A  - Transition probability matrix (K x K).
%   B  - Emission probability matrix (K x M).
%   pi - Initial state probability vector (K x 1).

fid = fopen(filename, 'r');
if fid == -1
    error('Cannot open file: %s', filename);
end

% Initialize variables.
states = {};
pi = [];
outputs = {};
transitions = {};

% Parse the file line by line.
while ~feof(fid)
    line = fgetl(fid);
    
    if ~ischar(line) || isempty(strtrim(line))
        continue;
    end
    
    % Split line into tokens.
    tokens = strsplit(strtrim(line));
    if isempty(tokens)
        continue;
    end
    
    keyword = tokens{1};
    
    if strcmp(keyword, 'State')
        % State line: State A 0.5
        stateName = tokens{2};
        initProb = str2double(tokens{3});
        states{end+1} = stateName; %#ok<AGROW>
        pi(end+1) = initProb; %#ok<AGROW>
        
    elseif strcmp(keyword, 'Output')
        % Output line: Output 0.31 0.48 0.21
        outputProbs = cellfun(@str2double, tokens(2:end));
        outputs{end+1} = outputProbs; %#ok<AGROW>
        
    elseif strcmp(keyword, 'Transition')
        % Transition line: Transition 0.9 0.1
        transProbs = cellfun(@str2double, tokens(2:end));
        transitions{end+1} = transProbs; %#ok<AGROW>
    end
end

fclose(fid);

% Convert cell arrays to matrices.
K = length(states);
M = length(outputs{1});

A = zeros(K, K);
B = zeros(K, M);

for i = 1:K
    A(i, :) = transitions{i};
    B(i, :) = outputs{i};
end

pi = pi(:);  % Convert to column vector.

end