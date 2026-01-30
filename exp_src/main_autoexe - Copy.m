% main_autoexe.m
% Automation script for HMM Baum-Welch training and Viterbi decoding.
% This script manages directory structures, executes C++ binaries, and aggregates results.

clear; clc; close all;

%=============================================================================-
% compile.

%% 1. Path Configuration.
% Get the directory where this current script is located (i.e., exp_src).
scriptPath = fileparts(mfilename('fullpath'));

% Navigate one level up to get the root directory.
rootDir = fileparts(scriptPath);

% Define the source and binary directories based on the root directory.
sourceDir = fullfile(rootDir, 'src');
binaryDir = fullfile(rootDir, 'bin');

% List the target file names (without extensions).
targetNames = {'BaumWelch', 'Viterbi'};

%% 2. Directory Preparation.
% Check if the binary directory exists.
if ~exist(binaryDir, 'dir')
    fprintf('Creating output directory: %s\n', binaryDir);
    mkdir(binaryDir);
end

%% 3. Compilation Loop.
% Loop through each target to compile them separately.
for i = 1:length(targetNames)
    baseName = targetNames{i};
    
    % Define the full path for the input C++ file.
    cppFile = fullfile(sourceDir, [baseName, '.cpp']);
    
    % Define the full path for the output Windows Executable file.
    exeFile = fullfile(binaryDir, [baseName, '.exe']);
    
    % Construct the compilation command for Windows (using g++).
    % Quotes are added to paths to handle potential spaces.
    compileCmd = sprintf('g++ -O3 "%s" -o "%s"', cppFile, exeFile);
    
    % Execute the command.
    fprintf('Compiling %s...\n', [baseName, '.cpp']);
    [status, cmdout] = system(compileCmd);
    
    % Check for compilation errors.
    if status ~= 0
        error('Compilation failed for %s. Output:\n%s', baseName, cmdout);
    else
        fprintf('Success: %s created.\n', [baseName, '.exe']);
        if ~isempty(cmdout)
            disp(cmdout);
        end
    end
end

disp('Build process finished.');
% =======================================================================

%% Configuration
% Define the target samples, state range, and number of seeds.
targetSamples = {'sample0'};
stateRange = 2:6;
numSeeds = 5; 
% Executable paths (Assuming compiled binaries are in bin folder).
exeBaumWelch = fullfile('..', 'bin', 'BaumWelch.exe');
exeViterbi = fullfile('..', 'bin', 'Viterbi.exe');

% Output directory base.
resultBaseDir = fullfile('..', 'results');

%% Main Loop for Training
% Iterate through each sample to train models.
for sIdx = 1:length(targetSamples)
    sampleName = targetSamples{sIdx};
    fprintf('Processing %s...\n', sampleName);

    % Define data directory path relative to execution folder.
    dataDir = fullfile('..', 'data', 'raw', 'samples', sampleName);
    
    % Create directory structure as defined in requirements.
    % results/[SampleName]/models
    % results/[SampleName]/history
    % results/[SampleName]/viterbi
    sampleDir = fullfile(resultBaseDir, sampleName);
    modelsDir = fullfile(sampleDir, 'models');
    historyDir = fullfile(sampleDir, 'history');
    viterbiDir = fullfile(sampleDir, 'viterbi');
    
    if ~exist(modelsDir, 'dir'), mkdir(modelsDir); end
    if ~exist(historyDir, 'dir'), mkdir(historyDir); end
    if ~exist(viterbiDir, 'dir'), mkdir(viterbiDir); end
    
    % Initialize master table data for this sample.
    % Columns: SampleName, States, Seed, TotalLen, FinalLogLik, Iterations, Params, AIC, BIC, ModelPath, HistoryPath
    masterData = {};
    
    for K = stateRange
        for seed = 1:numSeeds
            % Construct command for BaumWelch.exe.
            % Usage: [SampleName] [OutputDir] [States] [Seed]
            inputDir = fullfile('..', 'data', 'raw', 'samples', sampleName);
            cmd = sprintf('%s "%s" "%s" %d %d', exeBaumWelch, inputDir, modelsDir, K, seed);
            
            % Execute and capture stdout.
            [status, cmdout] = system(cmd);
            
            if status ~= 0
                warning('Execution failed for %s, K=%d, Seed=%d.', sampleName, K, seed);
                continue;
            end
            
            % Parse stdout for HISTORY and RESULTS.
            % History format: HISTORY [Step] [LogLik] [DError]
            historyLines = regexp(cmdout, 'HISTORY\s+(\d+)\s+([-\d\.]+)\s+([-\d\.]+)', 'tokens');
            
            % Results format: RESULTS [SampleName] [States] [Seed] [Symbols] [TotalLen] [FinalLogLik] [Iterations]
            resultLine = regexp(cmdout, 'RESULTS\s+(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+([-\d\.]+)\s+(\d+)', 'tokens');
            
            if ~isempty(resultLine)
                res = resultLine{1};
                % Parse numeric values.
                st = str2double(res{2});
                sd = str2double(res{3});
                M = str2double(res{4}); % Symbols
                T = str2double(res{5}); % TotalLength
                LL = str2double(res{6});
                Iter = str2double(res{7});
                
                % Calculate AIC and BIC.
                % Free parameters P = K*(K-1) + K*(M-1) + (K-1)
                % Transitions + Emissions + Initial
                P = st*(st-1) + st*(M-1) + (st-1);
                
                aic = -2 * LL + 2 * P;
                bic = -2 * LL + P * log(T);
                
                % Define file paths.
                modelFile = fullfile(modelsDir, sprintf('best_model_s%d_%d.txt', st, sd));
                histFile = fullfile(historyDir, sprintf('history_s%d_%d.tsv', st, sd));
                
                % Save history to TSV.
                fidHist = fopen(histFile, 'w');
                fprintf(fidHist, 'Step\tLogLik\tDError\n');
                for h = 1:length(historyLines)
                    hDat = historyLines{h};
                    fprintf(fidHist, '%s\t%s\t%s\n', hDat{1}, hDat{2}, hDat{3});
                end
                fclose(fidHist);
                
                % Append to master data.
                % SampleName, States, Seed, TotalLen, FinalLogLik, Iterations, Params, AIC, BIC, ModelPath, HistoryPath
                masterData(end+1, :) = {sampleName, st, sd, T, LL, Iter, P, aic, bic, modelFile, histFile};
            end
        end
    end
    
    % Append results to the global master.tsv.
    masterFile = fullfile(resultBaseDir, 'master.tsv');
    writeHeader = ~exist(masterFile, 'file');
    fidMaster = fopen(masterFile, 'a');
    if writeHeader
        fprintf(fidMaster, 'SampleName\tStates\tSeed\tTotalLen\tFinalLogLik\tIterations\tParams\tAIC\tBIC\tModelPath\tHistoryPath\n');
    end
    for r = 1:size(masterData, 1)
        fprintf(fidMaster, '%s\t%d\t%d\t%d\t%.4f\t%d\t%d\t%.4f\t%.4f\t%s\t%s\n', ...
            masterData{r,1}, masterData{r,2}, masterData{r,3}, masterData{r,4}, ...
            masterData{r,5}, masterData{r,6}, masterData{r,7}, masterData{r,8}, masterData{r,9}, ...
            masterData{r,10}, masterData{r,11});
    end
    fclose(fidMaster);
    if isempty(masterData)
        warning('No valid training data for %s. Skipping Viterbi.', sampleName);
        continue;
    end
    %% Viterbi Decoding using the Best Model
    % Find the model with Minimum BIC for this sample.
    sampleTable = cell2table(masterData, 'VariableNames', ...
        {'SampleName', 'States', 'Seed', 'TotalLen', 'FinalLogLik', 'Iterations', 'Params', 'AIC', 'BIC', 'ModelPath', 'HistoryPath'});
    
    [~, minIdx] = min(sampleTable.BIC);
    bestModel = sampleTable(minIdx, :);
    
    fprintf('Best Model for %s: States=%d, Seed=%d (BIC=%.2f). Running Viterbi...\n', ...
        sampleName, bestModel.States, bestModel.Seed, bestModel.BIC);
    
% Construct command for Viterbi.exe.
    % Usage: [SampleName] [OutputDir] [ModelPath] [States] [Seed]
    inputDir = fullfile('..', 'data', 'raw', 'samples', sampleName);
    
    cmdViterbi = sprintf('%s "%s" "%s" "%s" %d %d', ...
    exeViterbi, inputDir, viterbiDir, bestModel.ModelPath{1}, bestModel.States, bestModel.Seed);
    
    [vStatus, vOut] = system(cmdViterbi);
    if vStatus ~= 0
        warning('Viterbi execution failed for %s.', sampleName);
    else
        % We can optionally parse vOut if needed, but the C++ code handles file output.
        fprintf('Viterbi completed for %s.\n', sampleName);
    end
    
end

fprintf('All processes finished.\n');