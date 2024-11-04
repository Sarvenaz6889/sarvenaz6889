%% FM siusoid signal 3
%parameters
 %amplituse
A = 1;  
 %frequency
f0 = 5;         % Frequency in Hz
sigma = 0.1;    % Time constant (controls the Gaussian width)
phi = pi/4;     % Phase in radians
timeVec = 0:samplIntrvl:sigLen;
% Sine-Gaussian signal equation
sigVec =A *exp(-(timeVec.^2) / (2 * sigma^2)).*sin(2 *pi* f0 *timeVec + phi);

% Plot the signal
plot(timeVec, sigVec,'-');