%% Mock Data Challenge: Detect and Estimate a Quadratic Chirp Signal
clear all;
close all;

% Load the training data (noise-only) and the analysis data
load('/Users/sarvenazsadeghhasani/Downloads/TrainingData.mat'); % Noise-only data: trainData
load('/Users/sarvenazsadeghhasani/Downloads/analysisData.mat'); % Signal data: dataVec

% Sampling parameters
fs = sampFreq;             % Sampling frequency
nSamples = length(dataVec);% Number of samples
timeVec = (0:(nSamples-1)) / fs; % Time vector
SNR = 10;                     % Signal to noise ratio

% Estimate PSD using Welch's method
window = 256;               % Window size
noverlap = window / 2;      % 50% overlap
nfft = 2 * fs;              % Number of FFT points
[psd, freqVec] = pwelch(trainData, window, noverlap, nfft, fs);

% PSD Plot
figure; 
plot(freqVec, 10*log10(psd), 'LineWidth', 1.5);
xlabel('Frequency (Hz)'); ylabel('Power/Frequency (dB/Hz)');
title('Estimated PSD of Training Data (Noise)');
grid on;


psdPosFreq = psd'; % Interpolated PSD for positive frequencies
figure;
plot(timeVec, dataVec);
title('Analysis Data: Quadratic Chirp Signal');
xlabel('Time (s)'); ylabel('Amplitude');
grid on;

% Parameters for PSO optimization 
numRuns = 8;                  % Number of PSO runs
maxSteps = 2000;              % PSO iterations per run

searchRange = [40, 1, 1; ...  % Lower bounds: a1, a2, a3
               100, 50, 15];  % Upper bounds: a1, a2, a3

% Input Parameters for PSO
inParams = struct('dataX', timeVec, ...
                  'dataY', dataVec, ...
                  'dataXSq', timeVec.^2, ...
                  'dataXCb', timeVec.^3, ...
                  'snr', SNR, ...
                  'sampFreq', fs, ...
                  'psdPosFreq', psdPosFreq, ...
                  'rmin', searchRange(1, :), ...
                  'rmax', searchRange(2, :));
% Run PSO for Signal Detection 
disp('Starting PSO optimization for GLRT...');
rng('default');  % Set random seed for reproducibility
outStruct = pso4glrt(inParams, struct('maxSteps', maxSteps), numRuns);

% Plot Results
figure; 
hold on;
plot(timeVec, dataVec, 'k.', 'DisplayName', 'Input Data'); % Raw data

% Plot all runs' signals
for i = 1:length(outStruct.allRunsOutput)
    estSig = glrtSignal(outStruct.allRunsOutput(i).bestPos, timeVec);
    plot(timeVec, estSig, 'Color', [0.2, 0.8, 0.6], ...
        'DisplayName', ['Run ', num2str(i)]);
end

% Plot the best run's signal
plot(timeVec, outStruct.bestSig, 'r-', 'LineWidth', 2, ...
    'DisplayName', 'Best Fit Signal');
xlabel('Time (s)');
ylabel('Amplitude');
title('Estimated Signal Over PSO Runs');
legend('show');
grid on;

% Display Best Parameters
disp('Best Estimated Parameters for Quadratic Chirp:');
disp(['a1 = ', num2str(outStruct.bestQcCoefs(1))]);
disp(['a2 = ', num2str(outStruct.bestQcCoefs(2))]);
disp(['a3 = ', num2str(outStruct.bestQcCoefs(3))]);
