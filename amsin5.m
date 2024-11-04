%% AMSIUSOIDSIGNAL 5
%parameters

A = 1;           
f0 = 100;       
f1 = 5;      
phi = pi/4;      % Phase shift (radians)


sigLen= 1;
samplIntrvl= 0.008;
timeVec= 0:samplIntrvl:sigLen;


% Amplitude modulated signal
sigVec = A * cos(2 * pi * f1 * timeVec) .*(sin(f0 * timeVec+ phi));

% Plot the AM signal

plot(timeVec, sigVec, '-');
xlabel('Time (s)');
ylabel('Amplitude');
title('Amplitude Modulated (AM) Signal');
