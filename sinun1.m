%% 1 plot a sinusodial 
%signal parameters
A= 10;

% signal lenght
sigLen=1;
samplIntrvl=0.008;
timeVec= 0:samplIntrvl:sigLen;

% parameters for sinusodial signal
 %frequency in Hz
f0= 5; 
 %phase
phi=pi/4;

% time vector
sigVec= A * sin(2 * pi * f0 * timeVec + phi);

plot(timeVec,sigVec, '-');
title('sinusodial signal');