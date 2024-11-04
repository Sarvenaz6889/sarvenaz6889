%% AMSIUSOIDSIGNAL 5
%parameters
a1 = 5;
a2 = 2;
a3 = 2;
A = 10;
f0 = 50;           % Carrier frequency in Hz
t = 0:0.001:1;     % Time vector from 0 to 1 second with step size 0.00               
f1 = 5;  
phi = pi/4;      % Phase shift (radians)



sigLen= 1;
samplIntrvl= 0.008;
timeVec= 0:samplIntrvl:sigLen;
maxFreq = a1+2*a2+3*a3;
nyqFreq = 2*maxFreq;
samplFreq = 5*nyqFreq; 
samplIntrvl = 1/samplFreq;
nSamples = length(timeVec);


% Amplitude modulated signal
sigVec = A * cos(2 * pi * f1 * timeVec) .*(sin(f0 * timeVec+ phi));

% Plot the AM signal

plot(timeVec, sigVec, '-');
xlabel('Time (s)');
ylabel('Amplitude');
title('Amplitude Modulated (AM) Signal');
