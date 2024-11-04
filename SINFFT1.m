%% plot FM sinosuid 
% signal parameters
a1 = 5;
a2 = 2;
a3 = 2;
A = 10;
 %frequency in Hz
f0= 5; 
 %phase
phi=pi/4;
t = 0:0.001:1;     % Time vector from 0 to 1 second with step size 0.00


maxFreq = a1+2*a2+3*a3;
nyqFreq = 2*maxFreq;
samplFreq = 5*nyqFreq; 
samplIntrvl = 1/samplFreq;
sigLen= 1;
samplIntrvl= 0.008;
timeVec= 0:samplIntrvl:sigLen;
nSamples = length(timeVec);
% time vector
sigVec= A * sin(2 * pi * f0 * timeVec + phi);


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
xlabel('Frequency');
ylabel('|FFT|');
title('Periodogram');