%% AMFMSINUSOID% 6
%PARAMETERS
A = 1;         
f0 = 100;       
f1 = 5;        
b = 2;      
phi_m = pi/4;   % Phase of the modulating signal
a1 = 5;
a2 = 2;
a3 = 2;
t = 0:0.001:1;     % Time vector from 0 to 1 second with step size 0.00


sigLen= 1;
samplIntrvl= 0.008;
timeVec= 0:samplIntrvl:sigLen;

maxFreq = a1+2*a2+3*a3;
nyqFreq = 2*maxFreq;
samplFreq = 5*nyqFreq; 
samplIntrvl = 1/samplFreq;
sigLen= 1;
samplIntrvl= 0.008;
timeVec= 0:samplIntrvl:sigLen;
nSamples = length(timeVec);
% AM-FM Sinusoid
sigVec =A * cos(2*pi*f1*timeVec).*sin(2*pi*f0 *timeVec+b*cos(2*pi*f1*timeVec));

% Plot AM-FM signal
plot(timeVec, sigVec,'-');
xlabel('Time (s)');
ylabel('Amplitude');
title('AM-FM Sinusoid');
