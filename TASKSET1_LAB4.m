% parameters
snr = 10;
nSamples = 2048;
sampFreq = 1024;
timeVec = (0:(nSamples-1))/sampFreq;

% Generate the signal that is to be normalized
% Implementing Linear Transient Chirp signal directly
A = 1; % Amplitude value does not matter as it will be changed in the normalization
ta = 0.5; % Start time
td = 1.0; % Decay time constant
f0 = 50; % Starting frequency
f1 = 300; % Ending frequency
sigVec = zeros(size(timeVec));
for idx = 1:length(timeVec)
    t = timeVec(idx);
    if t >= ta
        f = f0 + (f1-f0)*(t-ta)/td;
        sigVec(idx) = A * exp(-(t-ta)/td) * sin(2*pi*f*(t-ta));
    end
end

% We will use the noise PSD used in colGaussNoiseDemo.m but add a constant
% to remove the parts that are zero.
noisePSD = @(f) (f>=100 & f<=300).*(f-100).*(300-f)/10000 + 1;

% Generate the PSD vector to be used in the normalization. Should be
% generated for all positive DFT frequencies. 
dataLen = nSamples/sampFreq;
kNyq = floor(nSamples/2)+1;
posFreq = (0:(kNyq-1))*(1/dataLen);
psdPosFreq = noisePSD(posFreq);
figure;
plot(posFreq,psdPosFreq);
axis([0,posFreq(end),0,max(psdPosFreq)]);
xlabel('Frequency (Hz)');
ylabel('PSD ((data unit)^2/Hz)');

% Calculation of the norm
% Implementing inner product calculation directly
fftSig = fft(sigVec);
dataLen = sampFreq*nSamples;
normSigSqrd = (1/dataLen)*sum(abs(fftSig(1:kNyq)).^2./psdPosFreq);

% Normalize signal to specified SNR
sigVec = snr*sigVec/sqrt(normSigSqrd);


% Obtain LLR values for multiple noise realizations
nH0Data = 1000;
llrH0 = zeros(1,nH0Data);
for lp = 1:nH0Data
    noiseVec = randn(size(sigVec)); % Simple white noise for demonstration
    fftNoise = fft(noiseVec);
    llrH0(lp) = (1/dataLen)*real(sum(conj(fftSig(1:kNyq)).*fftNoise(1:kNyq)./psdPosFreq));
end
%Obtain LLR for multiple data (=signal+noise) realizations
nH1Data = 1000;
llrH1 = zeros(1,nH1Data);
for lp = 1:nH0Data
    noiseVec = randn(size(sigVec)); % Simple white noise for demonstration
    dataVec = noiseVec + sigVec;
    fftData = fft(dataVec);
    llrH1(lp) = (1/dataLen)*real(sum(conj(fftSig(1:kNyq)).*fftData(1:kNyq)./psdPosFreq));
end

% Signal to noise ratio estimate
estSNR = (mean(llrH1)-mean(llrH0))/std(llrH0);

figure;
histogram(llrH0);
hold on;
histogram(llrH1);
xlabel('LLR');
ylabel('Counts');
title(['Estimated SNR = ',num2str(estSNR)]);

% A noise realization
figure;
plot(timeVec,noiseVec);
xlabel('Time (sec)');
ylabel('Noise');

% Plot the data (signal+noise) realization and the signal in the time domain
figure;
plot(timeVec,dataVec);
hold on;
plot(timeVec,sigVec);
xlabel('Time (sec)');
ylabel('Amplitude');
title('Data Realization and Signal in Time Domain');

% Plot the periodogram of the signal and data
figure;
[pxx, f] = periodogram(dataVec, [], [], sampFreq);
[pxxSig, fSig] = periodogram(sigVec, [], [], sampFreq);
plot(f, 10*log10(pxx));
hold on;
plot(fSig, 10*log10(pxxSig));
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
title('Periodogram of Signal and Data');

% Plot the spectrogram of the data
figure;
spectrogram(dataVec, hamming(256), 128, [], sampFreq, 'yaxis');
title('Spectrogram of Data');