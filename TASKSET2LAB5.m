% time vector dataX
dataX = (0:499) / 1024;  % Time vector from 0 to approximately 0.4883 seconds (500 samples)

%  parameter struct P 
P = struct('linear_coeff', 5.0, 'quadratic_coeff', 3.0, 'cubic_coeff', 1.0);

% Create foo_new
H = @(x) foo_new(dataX, x, P);

% DSNR values to plot
SNR_values = [10, 12, 15];

% for plotting
figure;
hold on;

% Loop through each SNR value and plot the signal time series
for i = 1:length(SNR_values)
    SNR = SNR_values(i);
    signal = H(SNR);  % Call the function handle with current SNR
    plot(dataX, signal, 'DisplayName', ['SNR = ' num2str(SNR)]);
end

% Add labels and title to the plot
xlabel('Time (seconds)');
ylabel('Amplitude');
title('Signal Time Series for Different SNR Values');
legend show;  % Show legend
grid on;      % Add grid to the plot
hold off;

% Function Definitions

function sigVec = foo_new(dataX, SNR, P)
    % Generate a quadratic chirp signal
    phaseVec = P.linear_coeff * dataX + P.quadratic_coeff * dataX.^2 + P.cubic_coeff * dataX.^3;
    sigVec = sin(2*pi*phaseVec);
    sigVec = SNR * sigVec / norm(sigVec);
end