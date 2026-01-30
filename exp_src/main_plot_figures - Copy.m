% main_plot_figures.m
% Generates Figures 1-4 as defined in figure.txt based on the results.
% Focuses on visualization for 'sample0' as requested.

clear; clc; close all;

%% Configuration
targetSample = 'sample0';
resultBaseDir = fullfile('..', 'results');
masterFile = fullfile(resultBaseDir, 'master.tsv');
figuresDir = fullfile(resultBaseDir, targetSample, 'figures');

if ~exist(figuresDir, 'dir'), mkdir(figuresDir); end

% Load Master Data.
opts = detectImportOptions(masterFile, 'FileType', 'text');
T = readtable(masterFile, opts);

% Filter for target sample.
T = T(strcmp(T.SampleName, targetSample), :);
states = unique(T.States);

%% Figure 1: Stability Analysis (Sorted Line Plot)
% X: Trial Index (Sorted), Y: Final Log-Likelihood.
figure('Name', 'Fig1_Stability', 'Color', 'w');
hold on;
colors = linspace(0.1, 0.8, length(states)); % Grayscale compatible

for i = 1:length(states)
    st = states(i);
    subT = T(T.States == st, :);
    ll = sort(subT.FinalLogLik, 'descend');
    
    plot(ll, 'LineWidth', 1.5, 'DisplayName', sprintf('K=%d', st), ...
        'Color', [colors(i) colors(i) colors(i)]);
end

title('Figure 1: Stability Analysis (Sorted Log-Likelihood)');
xlabel('Sorted Trial Index');
ylabel('Final Log-Likelihood');
legend('Location', 'SouthWest');
grid on;
hold off;
saveas(gcf, fullfile(figuresDir, 'Fig1_Stability.png'));

%% Figure 2: Learning Convergence
% X: Iterations, Y: Log Error (Log Scale).
figure('Name', 'Fig2_Convergence', 'Color', 'w');
hold on;

for i = 1:length(states)
    st = states(i);
    % Find best seed for this state (Max LogLik).
    subT = T(T.States == st, :);
    [~, idx] = max(subT.FinalLogLik);
    bestRow = subT(idx, :);
    
    % Load history file.
    histPath = bestRow.HistoryPath{1};
    if exist(histPath, 'file')
        histData = readtable(histPath, 'FileType', 'text');
        
        % DError can be negative or zero, handle for log plot.
        % The requirements ask for "Log Error". Usually this implies distance to convergence.
        % Here we plot the 'DError' value reported by C++ code.
        % We take abs because DError approaches 0.
        yData = abs(histData.DError);
        yData(yData == 0) = NaN; % Avoid log(0)
        
        semilogy(histData.Step, yData, 'LineWidth', 1.5, ...
            'DisplayName', sprintf('K=%d', st), ...
            'Color', [colors(i) colors(i) colors(i)]);
    end
end

title('Figure 2: Learning Convergence (DError)');
xlabel('Iterations');
ylabel('Log Error (abs(Likelihood Change))');
legend('Location', 'NorthEast');
grid on;
hold off;
saveas(gcf, fullfile(figuresDir, 'Fig2_Convergence.png'));

%% Figure 3: Model Selection
% X: States, Left Y: AIC/BIC, Right Y: Max Log-Likelihood.
figure('Name', 'Fig3_ModelSelection', 'Color', 'w');

bestAICs = zeros(length(states), 1);
bestBICs = zeros(length(states), 1);
maxLLs = zeros(length(states), 1);

for i = 1:length(states)
    st = states(i);
    subT = T(T.States == st, :);
    [maxLL, idx] = max(subT.FinalLogLik);
    maxLLs(i) = maxLL;
    bestAICs(i) = subT.AIC(idx);
    bestBICs(i) = subT.BIC(idx);
end

yyaxis left
plot(states, bestAICs, '-o', 'LineWidth', 2, 'DisplayName', 'AIC', 'Color', 'k');
hold on;
plot(states, bestBICs, '--s', 'LineWidth', 2, 'DisplayName', 'BIC', 'Color', [0.4 0.4 0.4]);
ylabel('Information Criteria (AIC / BIC)');
ax = gca; ax.YColor = 'k';

yyaxis right
plot(states, maxLLs, '-^', 'LineWidth', 1.5, 'DisplayName', 'LogLik', 'Color', [0.6 0.6 0.6]);
ylabel('Max Log-Likelihood');
ax = gca; ax.YColor = [0.4 0.4 0.4];

xlabel('Number of States (K)');
title('Figure 3: Model Selection');
legend('Location', 'Best');
grid on;
saveas(gcf, fullfile(figuresDir, 'Fig3_ModelSelection.png'));

%% Figure 4: Structure Analysis (Best Model)
% Select the global best model based on BIC.
[~, globalBestIdx] = min(T.BIC);
globalBest = T(globalBestIdx, :);

% Load HMM Parameters (A, B).
[A, B, Pi] = loadHMM(globalBest.ModelPath{1}, globalBest.States, 10); % Assuming 10 symbols max

% Load Viterbi Path and Observation.
% Path format: results/[SampleName]/viterbi/s[K]_[Seed]_out.txt
% Raw Observation: data/raw/samples/[SampleName]/1.txt (Assuming plotting first sample)
viterbiFile = fullfile(resultBaseDir, targetSample, 'viterbi', ...
    sprintf('s%d_%d_out.txt', globalBest.States, globalBest.Seed));
% Note: The path to raw data depends on the setup. Assuming standard path.
obsFile = fullfile('..', 'data', 'raw', 'samples', targetSample, '1.txt');

if exist(viterbiFile, 'file') && exist(obsFile, 'file')
    % Read files.
    % Viterbi output is characters 'A', 'B'... need conversion to 1, 2...
    vText = fileread(viterbiFile);
    vSeq = double(strsplit(strtrim(vText))) - double('A') + 1; % Convert A->1, B->2
    
    obsData = load(obsFile); % Assuming space delimited integers
    
    figure('Name', 'Fig4_Structure', 'Color', 'w', 'Position', [100 100 800 600]);
    
    % Subplot 1: Matrix A Heatmap.
    subplot(2, 2, 1);
    imagesc(A);
    colormap(gca, 'gray');
    colorbar;
    title('Transition Matrix A');
    xlabel('To State'); ylabel('From State');
    axis square;
    
    % Subplot 2: Matrix B Heatmap.
    subplot(2, 2, 2);
    imagesc(B);
    colormap(gca, 'gray');
    colorbar;
    title('Emission Matrix B');
    xlabel('Symbol'); ylabel('State');
    axis square;
    
    % Subplot 3: Viterbi Path (Time Series).
    subplot(2, 1, 2);
    hold on;
    % Plot observations (normalized to state range for overlay visualization).
    plot(obsData, 'k.', 'MarkerSize', 5, 'DisplayName', 'Observation');
    % Plot State Sequence.
    stairs(vSeq, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Viterbi State');
    
    xlabel('Time');
    ylabel('State / Symbol Index');
    title(sprintf('Viterbi Path (Best Model: K=%d)', globalBest.States));
    legend;
    ylim([0 max(max(obsData), globalBest.States)+1]);
    grid on;
    hold off;
    
    saveas(gcf, fullfile(figuresDir, 'Fig4_Structure.png'));
else
    warning('Could not find Viterbi output or Raw Observation for Figure 4.');
end

fprintf('Figure generation completed.\n');

%% Helper Function: Load HMM
function [A, B, Pi] = loadHMM(filepath, K, M)
    % Parses the specific HMM text format.
    fid = fopen(filepath, 'r');
    if fid == -1, error('Cannot open model file.'); end
    
    A = zeros(K, K);
    B = zeros(K, M); % Assuming max M symbols, will resize if needed
    Pi = zeros(K, 1);
    
    % Naive parsing loop based on format "State X Prob", "Output ...", "Transition ..."
    currentState = 0;
    while ~feof(fid)
        line = fgetl(fid);
        if isempty(line), continue; end
        
        if startsWith(line, 'State')
            % Parse: State [Char] [Prob]
            tokens = split(line);
            stateChar = tokens{2};
            stateIdx = double(stateChar) - double('A') + 1;
            Pi(stateIdx) = str2double(tokens{3});
            currentState = stateIdx;
        elseif startsWith(line, 'Output')
            % Parse: Output [p1] [p2] ...
            tokens = split(line);
            probs = str2double(tokens(2:end));
            B(currentState, 1:length(probs)) = probs;
        elseif startsWith(line, 'Transition')
            % Parse: Transition [p1] [p2] ...
            tokens = split(line);
            probs = str2double(tokens(2:end));
            A(currentState, 1:length(probs)) = probs;
        end
    end
    fclose(fid);
end