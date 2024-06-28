%%
clear all  % Clear all variables from the workspace.
close all  % Close all figure windows.
dir = "my working directory"
cd(dir)  % Change directory to the specified path.

%%
dirout = '..\Results_Simulations';  % Output directory for results.

fileresult = fullfile(dirout,strcat('Results_Correlated_ADD_MULT_', char(datetime("today")), '.mat'));
% Create the full file path for saving results, including today's date.

%% Population parameters

SIZE = [0.5 1 2 4 8 16];  % Different noise variances to be tested.

MaxRep = 250;  % Maximum number of repetitions for each simulation.

% Number of observations
nObs = 2^16;  % Number of observations in each simulated data set.

% Data distributions
DIST = [1 2 3 4 5];  % Different data distributions to be tested.

% Type of correlation matrix for the noise
NOISE_CORR_TYPE = [1 2 3 4];  % Different types of noise correlation matrices to be tested.

% Magnitude of correlation of the noise
NOISE_CORR_LEVEL = [2 4 6 8 10];  % Different levels of noise correlation magnitudes.

% Loading types
load('LOADING_TYPE_09Dec2023.mat')  % Load different loading types from file.

% Correlation matrices for noise
load('MATRIX_CORRELATIONS.mat')  % Load different noise correlation matrices from file.

% Total number of simulated data sets
omg = numel(DIST) * size(LOADING_TYPE,1) * numel(SIZE) * numel(NOISE_CORR_TYPE) * numel(NOISE_CORR_LEVEL) * MaxRep;
all_nobs = 0;  % Initialize the total count of observations.

% Number of variables and factors
nVar = 8;  % Number of variables in each data set.
Q = 3;  % Number of factors in the factor analysis.

% Simulations
for d = 1:numel(DIST)  % Loop over each data distribution.

    for corrtype = 1:size(LOADING_TYPE,1)  % Loop over each loading type.

        Lpop = LOADING_TYPE{corrtype,1};  % Load the population loading matrix.
        mylegend{corrtype} = LOADING_TYPE{corrtype,2};  % Load the legend for the loading type.

        for nct = 1:numel(NOISE_CORR_TYPE)  % Loop over each noise correlation type.

            if nct == 1
                ERRCORR = AVERAGE_CORR;  % Set the correlation matrix for average correlation.
            elseif nct == 2
                ERRCORR = HUB_CORR;  % Set the correlation matrix for hub correlation.
            elseif nct == 3
                ERRCORR = TOEPLITZ_CORR;  % Set the correlation matrix for Toeplitz correlation.
            elseif nct == 4
                ERRCORR = UNIFORM_CORR;  % Set the correlation matrix for uniform correlation.
            end

            for corr_level = 1:numel(NOISE_CORR_LEVEL)  % Loop over each noise correlation level.

                for i = 1:numel(SIZE)  % Loop over each noise variance.

                    ErrCorrMat = ERRCORR{corr_level,1};  % Select the appropriate error correlation matrix.

                    % Sample size
                    NoiseVariance = SIZE(i);  % Set the noise variance for this iteration.

                    ErrCorrSigma = corr2cov(sqrt(NoiseVariance * ones(1, nVar)), ErrCorrMat);
                    % Calculate the covariance matrix based on the noise variance and correlation matrix.

                    for r = 1:MaxRep  % Loop for repetitions.
                        all_nobs = all_nobs + 1;  % Increment the observation count.

                        disp(sprintf('%i %i %i %i %i %i %i %i', d, corrtype, i, r, nct, corr_level, all_nobs, omg))
                        % Display the current simulation parameters.

                        if d == 1
                            F0 = randn(nObs, nVar);  % Generate normally distributed data.
                        elseif d == 2
                            x = randn(nObs, nVar) * 0.8326 - 0.3466;
                            F0 = exp(x);  % Generate log-normal distributed data.
                        elseif d == 3
                            F = ramberg(0, -0.087, -0.0443, -0.0443, nObs, nVar);
                            F = F ./ sqrt(nObs);  % Generate Ramberg-Schmeiser distributed data.
                        elseif d == 4
                            F = ramberg(-0.886, 0.1333, 0.0193, 0.1588, nObs, nVar);
                            F = F ./ sqrt(nObs);  % Generate another type of Ramberg-Schmeiser distributed data.
                        elseif d == 5
                            a = -sqrt(12)/2; b = -a;
                            F = a + b * rand(nObs, nVar);  % Generate uniformly distributed data.
                        end

                        % Generate sample data noise-free
                        X0samp = F0 * Lpop';

                        % Generate Gaussian i.i.d. noise
                        NOISE = mvnrnd(zeros(1, nVar), ErrCorrSigma, nObs);

                        % Add noise (additive)
                        Xsamp_ADD = X0samp + NOISE;

                        % Calculate correlation matrix
                        Csamp_ADD = corr(Xsamp_ADD);

                        [~, Dsamp_ADD, Vsamp_ADD] = svd(Csamp_ADD);  % Singular Value Decomposition.

                        % The sample loadings are in Vsamp
                        Lsamp_ADD = Vsamp_ADD * sqrt(Dsamp_ADD);

                        Lsamp_ADD_ROT = rotatefactors(Lsamp_ADD(:, 1:Q), 'Method', 'procrustes', 'Target', Lpop(:, 1:Q), 'Type', 'Oblique');
                        % Rotate the factors using procrustes rotation.

                        %% Quality metrics
                        [PHI_ADD_ALL, ~] = congruence(Lsamp_ADD_ROT(:, 1:Q), Lpop(:, 1:Q));
                        % Calculate congruence between rotated sample loadings and population loadings.

                        for q = 1:Q
                            STAT_ADD = ReconstructionStatistics2(Lsamp_ADD_ROT(:, q), Lpop(:, q));
                            % Calculate reconstruction statistics for each component.

                            R2_PerComponent_ADD(1, q) = STAT_ADD.R2L;
                            MSE_PerComponent_ADD(1, q) = STAT_ADD.MSEL;
                            PHI_PerComponent_ADD(1, q) = PHI_ADD_ALL(q, 1);
                        end

                        PHI_ADD(all_nobs) = mean(abs(PHI_ADD_ALL));  % Mean congruence for all components.
                        MSE_ADD(all_nobs) = mean(MSE_PerComponent_ADD);  % Mean squared error for all components.

                        %%
                        % Groups for ANOVA
                        G_size(all_nobs, 1) = i;  % Group for noise variance.
                        G_corr(all_nobs, 1) = corrtype;  % Group for loading type.
                        G_err_type(all_nobs, 1) = nct;  % Group for noise correlation type.
                        G_corr_level(all_nobs, 1) = corr_level;  % Group for noise correlation level.
                        G_dist(all_nobs, 1) = d;  % Group for data distribution.

                        %% MULTIPLICATIVE NOISE

                        % Add noise (multiplicative)
                        Xsamp_MULT = X0samp + X0samp .* NOISE;

                        % Calculate correlation matrix
                        Csamp_MULT = corr(Xsamp_MULT);

                        [~, Dsamp_MULT, Vsamp_MULT] = svd(Csamp_MULT);  % Singular Value Decomposition.

                        % The sample loadings are in Vsamp
                        Lsamp_MULT = Vsamp_MULT * sqrt(Dsamp_MULT);

                        Lsamp_MULT_ROT = rotatefactors(Lsamp_MULT(:, 1:Q), 'Method', 'procrustes', 'Target', Lpop(:, 1:Q), 'Type', 'oblique');
                        % Rotate the factors using procrustes rotation.

                        %% Quality metrics
                        [PHI_MULT_ALL, ~] = congruence(Lsamp_MULT_ROT(:, 1:Q), Lpop(:, 1:Q));
                        % Calculate congruence between rotated sample loadings and population loadings.

                        for q = 1:Q
                            STAT_MULT = ReconstructionStatistics2(Lsamp_MULT_ROT(:, q), Lpop(:, q));
                            % Calculate reconstruction statistics for each component.

                            R2_PerComponent_MULT(1, q) = STAT_MULT.R2L;
                            MSE_PerComponent_MULT(1, q) = STAT_MULT.MSEL;
                            PHI_PerComponent_MULT(1, q) = PHI_MULT_ALL(q, 1);
                        end

                        PHI_MULT(all_nobs) = mean(abs(PHI_MULT_ALL));  % Mean congruence for all components.
                        MSE_MULT(all_nobs) = mean(MSE_PerComponent_MULT);  % Mean squared error for all components.

                    end
                end
            end  % Close noise correlation level loop.
        end  % Close noise correlation type loop.
    end  % Close loading type loop.
end  % Close data distribution loop.

%% Calculation ratio
% Calculate the ratio of MSE and PHI between additive and multiplicative noise.

FC_MSE = MSE_ADD ./ MSE_MULT;
FC_PHI = PHI_ADD ./ PHI_MULT;

RESULT.MSE_ADD = MSE_ADD;
RESULT.PHI_ADD = PHI_ADD;

RESULT.MSE_MULT = MSE_MULT;
RESULT.PHI_MULT = PHI_MULT;

RESULT.FC_MSE = FC_MSE;
RESULT.FC_PHI = FC_PHI;

% Create groups for ANOVA
GROUPS = {G_corr, G_dist, G_size, G_err_type, G_corr_level};
GROUP_names = {'Loading type', 'Distribution', 'Error variance', 'Error Corr struct', 'Error Corr level'};

RESULT.GROUPS = GROUPS;
RESULT.GROUP_names = GROUP_names;

% Save results
save(fileresult, 'RESULT', '-mat', '-v7.3')

disp('saved')  % Display 'saved' in Italian.
disp(strcat('fileresult =', '"', fileresult, '"'))  % Display the result file path.
