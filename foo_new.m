        % Function to generate a quadratic chirp signal
function sigVec = crcbgenqcsig(dataX, snr, P)
    % Generate a quadratic chirp signal
    % S = CRCBGENQSIG(X, SNR, P)
    % Generates a quadratic chirp signal S. 
    % X is the vector of time stamps at which the samples of the signal are to be computed. 
    % SNR is the matched filtering signal-to-noise ratio of S
    % P is a struct with fields:
    %   a1: coefficient for linear term
    %   a2: coefficient for quadratic term
    %   a3: coefficient for cubic term

    % Compute the phase of the signal using parameters from struct P
    phaseVec = P.linear_coeff * dataX + P.quadratic_coeff * dataX.^2 + P.cubic_coeff * dataX.^3;

    % Generate the signal
    sigVec = sin(2 * pi * phaseVec);

    % Normalize the signal to match the specified SNR
    sigVec = snr * sigVec / norm(sigVec);
end
function sigVec = foo_neww(dataX, SNR, P)
    % Generate a quadratic chirp signal
    % Inputs:
    %   dataX: Time samples vector
    %   SNR: Signal-to-noise ratio
    %   P: Struct containing signal parameters
    %     P.linear_coeff: Coefficient for linear term
    %     P.quadratic_coeff: Coefficient for quadratic term
    %     P.cubic_coeff: Coefficient for cubic term

    % Compute phase vector using struct fields
    phaseVec = P.linear_coeff * dataX + P.quadratic_coeff * dataX.^2 + P.cubic_coeff * dataX.^3;
    
    % Generate signal
    sigVec = sin(2*pi*phaseVec);
    
    % Normalize signal to specified SNR
    sigVec = SNR * sigVec / norm(sigVec);
end
% Example usage of the function
dataX = 0:0.01:10;  % Time vector from 0 to 10 seconds with increments of 0.01 seconds
snr = 10;           % Desired Signal-to-Noise Ratio

% Define parameters in a struct
P = struct('linear_coeff', 10, 'quadratic_coeff', 3, 'cubic_coeff', 3);  % Coefficients for the quadratic chirp

% Generate the signal using crcbgenqcsig
sigVec = crcbgenqcsig(dataX, snr, P);

% Plot the generated signal
figure;
plot(dataX, sigVec);
xlabel('Time (s)');
ylabel('Signal Amplitude');
title('Generated Quadratic Chirp Signal');
grid on;