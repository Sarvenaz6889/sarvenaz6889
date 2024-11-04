%% AMFMSINUSOID% 6
%PARAMETERS
A = 1;         
f0 = 100;       
f1 = 5;        
b = 2;      
phi_m = pi/4;   % Phase of the modulating signal

sigLen= 1;
samplIntrvl= 0.008;
timeVec= 0:samplIntrvl:sigLen;
% AM-FM Sinusoid
sigVec =A * cos(2*pi*f1*timeVec).*sin(2*pi*f0 *timeVec+b*cos(2*pi*f1*timeVec));

% Plot AM-FM signal
plot(timeVec, sigVec,'-');
xlabel('Time (s)');
ylabel('Amplitude');
title('AM-FM Sinusoid');
