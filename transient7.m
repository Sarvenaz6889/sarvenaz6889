
% Define parameters
A = 1;          % Amplitude of the chirp signal
t0 = 0;         % Start time of the chirp
f0 = 5;         % Initial frequency in Hz
f1 = 20;        % Final frequency in Hz
L = 1;          % Duration of the chirp in seconds

% Define the time vector
fs = 1000;      % Sampling frequency in Hz
sigLen= 1;
samplIntrvl= 0.008;
timeVec = 0:samplIntrvl:sigLen;  

% Chirp signal calculation for t in [t0, t0 + T]
sigvec = A * sin(2 * pi * (f0 * (timeVec) + f1 * (timeVec).^2));

% Plot the transient chirp signal
plot(timeVec, sigVec, '-');
