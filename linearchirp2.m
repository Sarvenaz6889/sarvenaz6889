%% linear chirp signal 2
% signal parameters
A=10; 

% signal length
sigLen= 1;
samplIntrvl= 0.008;
timeVec= 0:samplIntrvl:sigLen;

% signal and phase 
f0 = 5;         % Starting frequency in Hz
f1= 0.5;
k = 10;         % Chirp rate (rate of frequency change)
phi = pi/4;     % Phase in radians

 
sigVec= A * sin(2 * pi * (f0 * timeVec + f1 * k * timeVec.^2) + phi);


%plot the signal
plot(timeVec,sigVec, '-');


% plot the signal
figure;
plot(timeVec, sigVec,'-');
xlabel('Time (sec)');
title('Sampled signal');


%Length of data 
dataLen = timeVec(end)-timeVec(1);
%DFT sample corresponding to Nyquist frequency
kNyq = floor(nSamples/2)+1;
% Positive Fourier frequencies
posFreq = (0:(kNyq-1))*(1/dataLen);
% FFT of signal
fftSig = fft(sigVec);
% Discard negative frequencies
fftSig = fftSig(1:kNyq);

%Plot periodogram
figure;
    plot(posFreq,abs(fftSig));
title('Periodogram');