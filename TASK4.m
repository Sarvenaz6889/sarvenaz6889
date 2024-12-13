% Parameters
a1 = 10; a2 = 3; a3 = 3;
nSamples = 2048;
sampFreq = 1024;
timeVec = (0:(nSamples-1))' / sampFreq;  

% Noise PSD
noisePSD = @(f) (f >= 100 & f <= 300) .* (f - 100) .* (300 - f) / 10000 + 1;
dataLen = nSamples / sampFreq;
kNyq = floor(nSamples / 2) + 1;
posFreq = (0:(kNyq-1)) * (1 / dataLen);
psdPosFreq = noisePSD(posFreq);
data1 = randn(2048, 1); % Example data for data1.txt
data2 = randn(2048, 1); % Example data for data2.txt
data3 = randn(2048, 1); % Example data for data3.txt
save('DETEST/data1.txt', 'data1', '-ascii');
save('DETEST/data2.txt', 'data2', '-ascii');
save('DETEST/data3.txt', 'data3', '-ascii');

% Calculate GLRT for each data realization
glrtValues = zeros(1, 3);
for n = 1:3
    try
        dataVec = load(['DETEST/data', num2str(n), '.txt'], '-ascii');
        dataVec = dataVec(:);
        glrtValues(n) = glrtqcsig(dataVec, timeVec, psdPosFreq, [a1, a2, a3]);

    catch ME
        fprintf('Error processing data%d.txt: %s\n', n, ME.message);
    end
end

% M noise realizations under H₀ and calculate GLRT values
M = 10000;
glrtH0 = zeros(1, M);
for m = 1:M
    try
        noiseVec = statgaussnoisegen(nSamples, [posFreq(:), psdPosFreq(:)], m, sampFreq);
        glrtH0(m) = glrtqcsig(noiseVec, timeVec, psdPosFreq, [a1, a2, a3]);
    catch ME
        fprintf('Error in iteration %d: %s\n', m, ME.message);
    end
end

% Estimate significance for each data realization
significance = zeros(1, 3);
for n = 1:3
    significance(n) = sum(glrtH0 >= glrtValues(n)) / M;
end


% Display results
for n = 1:3
    fprintf('Data Realization %d:\n', n);
    fprintf('  GLRT Value: %.4f\n', glrtValues(n));
    fprintf('  Estimated Significance: %.6f\n\n', significance(n));
end

function sigVec = crcbgenqcsig(dataX, snr, qcCoefs)
   
    phaseVec = qcCoefs(1)*dataX + qcCoefs(2)*dataX.^2 + qcCoefs(3)*dataX.^3;
    sigVec = sin(2*pi*phaseVec)
    fprintf('Generated signal (sigVec) - first 10 values: \n');
    disp(sigVec(1:10));  % Display first 10 values of the signal
    sigVec = snr * sigVec / norm(sigVec)
    fprintf('Normalized signal (sigVec) - first 10 values: \n');
    disp(sigVec(1:10));  % Display first 10 values of the normalized signal
end

function [normSig, normFac] = normsig4psd(sig, sampFreq, psdVals, snr)
    fftSig = fft(sig);
    kNyq = floor(length(sig) / 2) + 1;
    fftSig = fftSig(1:kNyq);
    sigPower = sum(abs(fftSig).^2 ./ psdVals);
    normFac = snr / sqrt(sigPower * sampFreq / length(sig));
    normSig = normFac * sig;
end

% Inner product
function innProd = innerprodpsd(x, y, sampFreq, psdVals)
    fftX = fft(x);
    fftY = fft(y);
    kNyq = floor(length(x) / 2) + 1;
    
    % Debugging: Print the first few FFT values of x and y
    %fprintf('FFT of x (first few):\n');
    %disp(fftX(1:10));  % Display first 10 FFT values of x
    %fprintf('FFT of y (first few):\n');
    %disp(fftY(1:10));  % Display first 10 FFT values of y
    
    innProd = sampFreq * sum((fftX(1:kNyq) .* conj(fftY(1:kNyq))) ./ psdVals) / length(x);
end

function noiseVec = statgaussnoisegen(nSamples, psdVals, seed, sampFreq)
    rng(seed); % Set random seed
    noiseVec = randn(1, nSamples); % Generate white noise
    fftNoise = fft(noiseVec);
    
    % Interpolate PSD to match frequency resolution
    kNyq = floor(nSamples / 2) + 1;
    posFreq = (0:kNyq-1) * (sampFreq / nSamples);
    interpPSD = interp1(psdVals(:,1), psdVals(:,2), posFreq, 'linear', 'extrap');

    fftNoise(1:kNyq) = fftNoise(1:kNyq) .* sqrt(interpPSD * sampFreq / nSamples);
    fftNoise(kNyq+1:end) = conj(fftNoise(kNyq-1:-1:2)); % Hermitian symmetry
    noiseVec = real(ifft(fftNoise)); % Return time-domain noise
    noiseVec = noiseVec(:);
    function noiseVec = statgaussnoisegen(nSamples, psdVals, seed, sampFreq)
   
end

end
    % Debugging: Check the value of llr
    %fprintf('LLR Value: %.4f\n', llr);
    
function glrtValue = glrtqcsig(dataVec, timeVec, psdPosFreq, qcCoefs)
    dataVec = dataVec(:);
    timeVec = timeVec(:);
    psdPosFreq = psdPosFreq(:);
    
    sigVec = crcbgenqcsig(timeVec, 1, qcCoefs);
    sampFreq = 1 / (timeVec(2) - timeVec(1));
    [templateVec, ~] = normsig4psd(sigVec, sampFreq, psdPosFreq, 1);
    llr = innerprodpsd(dataVec, templateVec, sampFreq, psdPosFreq);
    glrtValue = llr^2;
end
  
% Estimate significance with increasing M
M_values = [1000, 5000, 10000, 20000];
significance = zeros(3, length(M_values));

for idx = 1:length(M_values)
    M = M_values(idx);
    glrtH0 = zeros(1, M);
    
    % Generate noise realizations
    for m = 1:M
        try
            noiseVec = statgaussnoisegen(nSamples, [posFreq(:), psdPosFreq(:)], m, sampFreq);
            glrtH0(m) = glrtqcsig(noiseVec, timeVec, psdPosFreq, [a1, a2, a3]);
        catch ME
            fprintf('Error in iteration %d: %s\n', m, ME.message);
        end
    end
    
    % Calculate significance for each realization
    for n = 1:3
        significance(n, idx) = sum(glrtH0 >= glrtValues(n)) / M;
    end
    
    fprintf('M = %d: Significance values = [%.6f %.6f %.6f]\n', ...
        M, significance(:, idx));
end

profile on
M = 10000;
glrtH0 = zeros(1, M);
for m = 1:M
    try
        noiseVec = statgaussnoisegen(nSamples, [posFreq(:), psdPosFreq(:)], m, sampFreq);
        glrtH0(m) = glrtqcsig(noiseVec, timeVec, psdPosFreq, [a1, a2, a3]);
    catch ME
        fprintf('Error in iteration %d: %s\n', m, ME.message);
    end
end
profile viewer
profile off