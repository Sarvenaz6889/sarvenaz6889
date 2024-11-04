%% FM sinosuid 4
% parameters
A = 1;             % Amplitude
f0 = 50;           % Carrier frequency in Hz
fm = 5;            % Modulation frequency in Hz
beta = 10;         % Frequency deviation (modulation index)
t = 0:0.001:1;     % Time vector from 0 to 1 second with step size 0.00

sigLen= 1;
samplIntrvl= 0.008;
timeVec= 0:samplIntrvl:sigLen;

% FM sinusoid signal equation
sigVec = A * sin(2 * pi * f0 * timeVec + beta * cos(2 * pi * fm * timeVec));

% Plot the FM sinusoid signal
figure;
plot(timeVec, sigVec,'-');
title('Frequency Modulated Sinusoid');
