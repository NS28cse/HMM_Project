% main_autoexe.m
% Automation script for HMM Baum-Welch training, model selection,
% and Viterbi decoding.

clear; clc; close all;

%% Experiment configuration.
samples   = {'sample0', 'sample1', 'sample2', 'sample3', 'sample4'};
stateList = 2:6;
% Seed specification.
% Examples:
%   seedList = 1:50;        % range
%   seedList = 1:2:49;      % step
%   seedList = [3 7 42];    % explicit list
%   seedList = 17;          % single seed
seedList = 1:50;

%% Path configuration.
scriptDir = fileparts(mfilename('fullpath'));
rootDir   = fileparts(scriptDir);

srcDir     = fullfile(rootDir, 'src');
binDir     = fullfile(rootDir, 'bin');
dataDir    = fullfile(rootDir, 'data', 'raw', 'samples');
resultsDir = fullfile(rootDir, 'results');

exeBaumWelch = fullfile(binDir, 'BaumWelch.exe');
exeViterbi   = fullfile(binDir, 'Viterbi.exe');

%% Compile C++ sources.
targets = {'BaumWelch', 'Viterbi'};

if ~exist(binDir, 'dir')
    mkdir(binDir);
end

for i = 1:numel(targets)
    cppFile = fullfile(srcDir, [targets{i}, '.cpp']);
    exeFile = fullfile(binDir, [targets{i}, '.exe']);

    cmd = sprintf('g++ -O3 "%s" -o "%s"', cppFile, exeFile);
    fprintf('Compiling %s.\n', cppFile);

    [status, out] = system(cmd);
    if status ~= 0
        error('Compilation failed:\n%s', out);
    end
end

%% Master TSV initialization.
masterPath = fullfile(resultsDir, 'master.tsv');
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end

fidMaster = fopen(masterPath, 'w');
fprintf(fidMaster, ...
    'SampleName\tStates\tSeed\tSymbols\tTotalLen\tFinalLogLik\tIterations\tParams\tAIC\tBIC\tModelPath\tHistoryPath\n');
fclose(fidMaster);

%% Training loop.
for s = 1:numel(samples)
    sampleName = samples{s};
    fprintf('Processing sample: %s\n', sampleName);

    sampleDataDir = fullfile(dataDir, sampleName);
    sampleResDir  = fullfile(resultsDir, sampleName);

    modelDir   = fullfile(sampleResDir, 'models');
    histDir    = fullfile(sampleResDir, 'history');
    viterbiDir = fullfile(sampleResDir, 'viterbi');

    if ~exist(modelDir, 'dir'), mkdir(modelDir); end
    if ~exist(histDir, 'dir'), mkdir(histDir); end
    if ~exist(viterbiDir, 'dir'), mkdir(viterbiDir); end

    for K = stateList
        bestLogLik = -inf;
        bestSeed   = NaN;
        bestModel  = '';

        for seed = seedList
            fprintf('Running BaumWelch: %s States=%d Seed=%d\n',sampleName, K, seed);

            cmd = sprintf('"%s" "%s" "%s" %d %d', ...
                exeBaumWelch, sampleDataDir, modelDir, K, seed);

            [status, out] = system(cmd);
            if status ~= 0
                warning('BaumWelch failed: K=%d Seed=%d.', K, seed);
                continue;
            end

            %% Save HISTORY.
            histFile = fullfile(histDir, sprintf('history_s%d_%d.tsv', K, seed));
            fidHist = fopen(histFile, 'w');
            fprintf(fidHist, 'Step\tLogLik\tDError\n');

            lines = strsplit(out, '\n');

            for iLine = 1:numel(lines)
                line = strtrim(lines{iLine});
              if strncmp(line, 'HISTORY', 7)
                vals = strsplit(line);
                fprintf(fidHist, '%s\t%s\t%s\n', vals{2}, vals{3}, vals{4});
              end
            end

            fclose(fidHist);

            %% Parse RESULTS.
            for iLine = 1:numel(lines)
              line = strtrim(lines{iLine});
              if strncmp(line, 'RESULTS', 7)
              vals = strsplit(line);

                Symbols     = str2double(vals{5});
                TotalLen    = str2double(vals{6});
                FinalLogLik = str2double(vals{7});
                Iterations  = str2double(vals{8});


                    fprintf('>RESULTS: LogLik=%.6g, Iter=%d\n', FinalLogLik, Iterations);


                    % Parameter count: pi + A + B.
                    params = (K - 1) + K * (K - 1) + K * (Symbols - 1);
                    AIC = -2 * FinalLogLik + 2 * params;
                    BIC = -2 * FinalLogLik + params * log(TotalLen);

                    modelPath = fullfile(modelDir, sprintf('markov_output_s%d_%d.txt', K, seed));

                    fidMaster = fopen(masterPath, 'a');
                    fprintf(fidMaster, ...
                        '%s\t%d\t%d\t%d\t%d\t%.6g\t%d\t%d\t%.6g\t%.6g\t%s\t%s\n', ...
                        sampleName, K, seed, Symbols, TotalLen, ...
                        FinalLogLik, Iterations, params, AIC, BIC, ...
                        modelPath, histFile);
                    fclose(fidMaster);

                    if FinalLogLik > bestLogLik
                        bestLogLik = FinalLogLik;
                        bestSeed   = seed;
                        bestModel  = modelPath;
                    end
                end
            end
        end

        %% Save best model.
        if ~isnan(bestSeed)
            bestModelPath = fullfile(modelDir, sprintf('best_model_s%d.txt', K));
            copyfile(bestModel, bestModelPath);

            %% Run Viterbi using best model.
            fprintf('Running Viterbi: %s States=%d BestSeed=%d\n',sampleName, K, bestSeed);

            cmd = sprintf('"%s" "%s" "%s" "%s" %d %d', ...
                exeViterbi, sampleDataDir, viterbiDir, bestModelPath, K, bestSeed);
            [status, out] = system(cmd);
        end
        files = dir(fullfile(viterbiDir, '*_out.txt'));
        fprintf('>Viterbi finished (%d sequences)\n', numel(files))
    end
end

fprintf('All processes finished successfully.\n');
