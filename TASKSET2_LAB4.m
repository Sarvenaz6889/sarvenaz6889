    % Parameters
snr = 10;
nSamples = 2048;
sampFreq = 1024;
timeVec = (0:(nSamples-1))/sampFreq;
% Generate the Linear Transient Chirp signal
A = 1; 
ta = 0.5;
td = 1.0;
f0 = 50;
f1 = 300;
sigVec = zeros(size(timeVec));
for idx = 1:length(timeVec)
    t = timeVec(idx);
    if t >= ta
        f = f0 + (f1-f0)*(t-ta)/td;
        sigVec(idx) = A * exp(-(t-ta)/td) * sin(2*pi*f*(t-ta));
    end
end

% Interpolate initial LIGO PSD
iLIGOdata = load('iLIGOSensitivity.txt');
freq = iLIGOdata(:,1);
psd = iLIGOdata(:,2);

dataLen = nSamples/sampFreq;
kNyq = floor(nSamples/2)+1;
posFreq = (0:(kNyq-1))*(1/dataLen);

psdPosFreq = interp1(freq, psd, posFreq, 'linear', 'extrap');

figure;
loglog(posFreq, psdPosFreq);
xlabel('Frequency (Hz)');
ylabel('PSD (1/Hz)');
title('Interpolated initial LIGO PSD');

% Calculation of the norm
fftSig = fft(sigVec);
dataLen = sampFreq*nSamples;
normSigSqrd = (1/dataLen)*sum(abs(fftSig(1:kNyq)).^2./psdPosFreq);
sigVec = snr*sigVec/sqrt(normSigSqrd);

% Test
nH0Data = 1000;
llrH0 = zeros(1,nH0Data);
llrH1 = zeros(1,nH0Data);

for lp = 1:nH0Data
    noiseVec = randn(size(sigVec)) .* sqrt(psdPosFreq(1));  % Simplified noise generation
    fftNoise = fft(noiseVec);
    llrH0(lp) = (1/dataLen)*real(sum(conj(fftSig(1:kNyq)).*fftNoise(1:kNyq)./psdPosFreq));
    
    dataVec = noiseVec + sigVec;
    fftData = fft(dataVec);
    llrH1(lp) = (1/dataLen)*real(sum(conj(fftSig(1:kNyq)).*fftData(1:kNyq)./psdPosFreq));
end

estSNR = (mean(llrH1)-mean(llrH0))/std(llrH0);

figure;
histogram(llrH0);
hold on;
histogram(llrH1);
xlabel('LLR');
ylabel('Counts');
title(['Estimated SNR = ',num2str(estSNR)]);

% Plot realizations
figure;
plot(timeVec,noiseVec);
xlabel('Time (sec)');
ylabel('Noise');
title('Noise Realization');

figure;
plot(timeVec,dataVec);
hold on;
plot(timeVec,sigVec);
xlabel('Time (sec)');
ylabel('Amplitude');
title('Data Realization and Signal');

disp(['Target SNR: ', num2str(snr)]);
disp(['Estimated SNR: ', num2str(estSNR)]);