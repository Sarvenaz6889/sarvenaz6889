%% linear chirp signal 2
% signal parameters
A=10; 
a1 = 5;
a2 = 2;
a3 = 2;
t = 0:0.001:1;     % Time vector from 0 to 1 second with step size 0.00
f0 = 5;         % Starting frequency in Hz
f1= 0.5;
k = 10;         % Chirp rate (rate of frequency change)
phi = pi/4;     % Phase in radians

% signal length
sigLen= 1;
samplIntrvl= 0.008;
timeVec= 0:samplIntrvl:sigLen;


sigVec= A * sin(2 * pi * (f0 * timeVec + f1 * k * timeVec.^2) + phi);


%plot the signal
plot(timeVec,sigVec, '-');

maxFreq = a1+2*a2+3*a3;
nyqFreq = 2*maxFreq;
samplFreq = 5*nyqFreq; 
samplIntrvl = 1/samplFreq;
sigLen= 1;
samplIntrvl= 0.008;
timeVec= 0:samplIntrvl:sigLen;
nSamples = length(timeVec);
% FM sinusoid signal equation
sigVec = A * sin(2 * pi * f0 * timeVec + beta * cos(2 * pi * fm * timeVec));
