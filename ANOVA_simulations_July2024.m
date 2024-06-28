% Change the current directory
cd('...');

% Clear all variables from the workspace
clear all;

% Define directories for simulations and output
dirsim = '...\Results_Simulations'; %directory where results of simulations are
dirout = '...\Results_ANOVA';

% Define the paths to the result files
fileresult1 = fullfile(dirsim, "Results_Uncorrelated_ADD_MULT_05-Jan-2024.mat"); %files from the simulations
fileresult2 = fullfile(dirsim, "Results_Correlated_ADD_MULT_05-Jan-2024.mat"); %files from the simulations

% Select responses to analyze
response1 = 'MSE_ADD';
response2 = 'MSE_MULT';

% Uncomment to analyze different responses
% response1 = 'PHI_ADD';
% response2 = 'PHI_MULT';

% response1 = 'FC_MSE';
% response2 = 'FC_MULT';

% Load results for uncorrelated noise simulations
load(fileresult1); 
R = RESULT;

% Extract the specified responses
y1 = R.(response1);
y2 = R.(response2);

% Perform ANOVA for uncorrelated additive noise simulations
% Rank transform the data
[~, anovaTableAddUncorr, anovaStatsAddUncorr] = anovan(rank_transform(y1), cell2mat(R.GROUPS), 'model', 'interaction', 'varnames', R.GROUP_names, 'display', 'off');

% Do not transform the data (commented out)
% [~, anovaTableAddUncorr, anovaStatsAddUncorr] = anovan(y1, cell2mat(R.GROUPS), 'model', 'interaction', 'varnames', R.GROUP_names, 'display', 'off');

% Calculate Eta squared for the ANOVA table
[anovaTableEtaAddUncorr, ~] = GetEtaSquare(anovaTableAddUncorr, 1);

% Extract date from the filename
fileresult1 = char(fileresult1);
date1 = fileresult1(end-14:end-4);

% Save ANOVA results for uncorrelated additive noise
file1 = fullfile(dirout, strcat(sprintf('ANOVA_%s_Uncorrelated_', response1), date1, '.xlsx'));
xlswrite(file1, table2cell(anovaTableEtaAddUncorr));
disp(sprintf('Written: %s', file1));
close all;

% Save ANOVA statistics
file1anova = fullfile(dirout, strcat(sprintf('ANOVA_STATS_%s_Uncorrelated_', response1), date1, '.mat'));
save(file1anova, "anovaStatsAddUncorr", '-mat');
disp(sprintf('Written: %s', file1anova));

% Save ANOVA table
file1anova2 = fullfile(dirout, strcat(sprintf('ANOVA_TABLE_%s_Uncorrelated_', response1), date1, '.mat'));
save(file1anova2, "anovaTableAddUncorr", '-mat');
disp(sprintf('Written: %s', file1anova2));


% Perform ANOVA for uncorrelated multiplicative noise simulations
[~, anovaTableMultUncorr, anovaStatsMultUncorr] = anovan(rank_transform(y2), cell2mat(R.GROUPS), 'model', 'interaction', 'varnames', R.GROUP_names, 'display', 'off');
[anovaTableEtaMultUncorr, ~] = GetEtaSquare(anovaTableMultUncorr, 1);

% Save ANOVA results for uncorrelated multiplicative noise
file2 = fullfile(dirout, strcat(sprintf('ANOVA_%s_Uncorrelated_', response2), date1, '.xlsx'));
xlswrite(file2, table2cell(anovaTableEtaMultUncorr));
disp(sprintf('Written: %s', file2));
close all;

% Save ANOVA statistics
file2anova = fullfile(dirout, strcat(sprintf('ANOVA_STATS_%s_Uncorrelated_', response2), date1, '.mat'));
save(file2anova, "anovaStatsMultUncorr", '-mat');
disp(sprintf('Written: %s', file2anova));

% Save ANOVA table
file2anova2 = fullfile(dirout, strcat(sprintf('ANOVA_TABLE_%s_Uncorrelated_', response2), date1, '.mat'));
save(file2anova2, "anovaTableMultUncorr", '-mat');
disp(sprintf('Written: %s', file2anova2));



% Clear variables for next load
clear R RESULTS;

% Load results for correlated noise simulations
load(fileresult2); 
R = RESULT;

% Extract the specified responses
y3 = R.(response1); % Additive
y4 = R.(response2); % Multiplicative

% Perform ANOVA for correlated additive noise simulations
[~, anovaTableAddCorr, anovaStatsAddCorr] = anovan(rank_transform(y3), cell2mat(R.GROUPS), 'model', 'interaction', 'varnames', R.GROUP_names, 'display', 'off');
[anovaTableEtaAddCorr, ~] = GetEtaSquare(anovaTableAddCorr, 1);

% Extract date from the filename
fileresult2 = char(fileresult2);
date2 = fileresult2(end-14:end-4);

% Save ANOVA results for correlated additive noise
file3 = fullfile(dirout, strcat(sprintf('ANOVA_%s_Correlated_', response1), date2, '.xlsx'));
xlswrite(file3, table2cell(anovaTableEtaAddCorr));
disp(sprintf('Written: %s', file3));
close all;

% Save ANOVA statistics
file3anova = fullfile(dirout, strcat(sprintf('ANOVA_STATS_%s_Correlated_', response1), date2, '.mat'));
save(file3anova, "anovaStatsAddCorr", '-mat');
disp(sprintf('Written: %s', file3anova));

% Save ANOVA table
file3anova2 = fullfile(dirout, strcat(sprintf('ANOVA_TABLE_%s_Correlated_', response1), date2, '.mat'));
save(file3anova2, "anovaTableAddCorr", '-mat');
disp(sprintf('Written: %s', file3anova2));


% Perform ANOVA for correlated multiplicative noise simulations
[~, anovaTableMultCorr, anovaStatsMultCorr] = anovan(rank_transform(y4), cell2mat(R.GROUPS), 'model', 'interaction', 'varnames', R.GROUP_names, 'display', 'off');
[anovaTableEtaMultCorr, ~] = GetEtaSquare(anovaTableMultCorr, 1);

% Save ANOVA results for correlated multiplicative noise
file4 = fullfile(dirout, strcat(sprintf('ANOVA_%s_Correlated_', response2), date2, '.xlsx'));
xlswrite(file4, table2cell(anovaTableEtaMultCorr));
disp(sprintf('Written: %s', file4));
close all;

% Save ANOVA statistics
file4anova = fullfile(dirout, strcat(sprintf('ANOVA_STATS_%s_Correlated_', response2), date2, '.mat'));
save(file4anova, "anovaStatsMultCorr", '-mat');
disp(sprintf('Written: %s', file4anova));

% Save ANOVA table
file4anova2 = fullfile(dirout, strcat(sprintf('ANOVA_TABLE_%s_Correlated_', response2), date2, '.mat'));
save(file4anova2, "anovaTableMultCorr", '-mat');
disp(sprintf('Written: %s', file4anova2));


