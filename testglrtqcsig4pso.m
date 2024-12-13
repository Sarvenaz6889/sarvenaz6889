sampFreq = 1024; % Hz
sampTime = 5; % seconds
timeVec = 0:1/sampFreq:sampTime-1/sampFreq;

a1_0 = 10; 
a2_0 = 3;
a3_0 = 0.5;
qcSig = crcbgenqcsig(timeVec, 1, [a1_0, a2_0, a3_0]);

% Generate colored noise 
psdFunc = @(f) (f >= 0 & f <= sampFreq/2).*(1./(f + 1));
psdPosFreq = psdFunc(0:sampFreq/length(timeVec):sampFreq/2);
% Normalize signal to SNR = 10
[qcSig, ~] = normsig4psd(qcSig, sampFreq, psdPosFreq, 10);
% Generate colored noise
noiseVec = statgaussnoisegen(length(timeVec), sampFreq, psdPosFreq, 10);
% Create data realization
dataY = qcSig + noiseVec;
% array of values for a1
a1_min = 5;
a1_max = 15;
delta_a1 = 0.1;
A = a1_min:delta_a1:a1_max;

% Ranges for a2 and a3
a2_min = 1;
a2_max = 5;
a3_min = 0;
a3_max = 1;

% Standardized coordinate values
X = zeros(length(A), 3);
X(:,1) = (A - a1_min) / (a1_max - a1_min);
X(:,2) = (a2_0 - a2_min) / (a2_max - a2_min);
X(:,3) = (a3_0 - a3_min) / (a3_max - a3_min);

% Fitness values using glrtqcsig4pso
P.dataY = dataY;
P.dataX = timeVec;
P.dataXSq = timeVec.^2;
P.dataXCb = timeVec.^3;
P.psdPosFreq = psdPosFreq;
P.rmin = [a1_min, a2_min, a3_min];
P.rmax = [a1_max, a2_max, a3_max];

[fitnessVals, ~] = glrtqcsig4pso(X, P);

% Plot fitness values against A
figure;
plot(A, fitnessVals);
xlabel('a1');
ylabel('Fitness Value');
title('Fitness Values vs a1');
grid on;

% min
[minFitness, minIdx] = min(fitnessVals);
hold on;
plot(A(minIdx), minFitness, 'ro', 'MarkerSize', 10);
text(A(minIdx), minFitness, sprintf('  Min at a1 ≈ %.2f', A(minIdx)), 'VerticalAlignment', 'bottom');

% Add a1 
line([a1_0, a1_0], ylim, 'Color', 'r', 'LineStyle', '--');
text(a1_0, max(ylim), sprintf('True a1 = %.2f', a1_0), 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');

legend('Fitness Values', 'Global Minimum', 'True a1');

% Function definitions
function qcSig = crcbgenqcsig(timeVec, amp, params)
    phaseVec = params(1)*timeVec + params(2)*timeVec.^2 + params(3)*timeVec.^3;
    qcSig = amp * sin(2*pi*phaseVec);
end

function [normSigVec, normFac] = normsig4psd(sigVec, sampFreq, psdVec, snr)
    nSamples = length(sigVec);
    kNyq = floor(nSamples/2)+1;
    normSigSqrd = innerprodpsd(sigVec, sigVec, sampFreq, psdVec);
    normFac = snr/sqrt(normSigSqrd);
    normSigVec = normFac*sigVec;
end

function noiseVec = statgaussnoisegen(nSamples, sampFreq, psdVec, seed)
    rng(seed);
    noiseVec = randn(1, nSamples);
    noiseVec = colorGaussianNoise(noiseVec, sampFreq, psdVec);
end

function coloredNoise = colorGaussianNoise(whiteNoise, sampFreq, psdVec)
    N = length(whiteNoise);
    df = sampFreq/N;
    fVec = (0:N-1)*df;
    sqrtPSD = sqrt(interp1(0:df:(sampFreq/2), psdVec, fVec, 'linear', 'extrap'));
    coloredNoise = real(ifft(fft(whiteNoise) .* sqrtPSD));
end

function [F, R] = glrtqcsig4pso(X, P)
    [nRows, ~] = size(X);
    F = zeros(nRows, 1);
    R = zeros(nRows, 3);
    for lpc = 1:nRows
        x = X(lpc,:).*(P.rmax - P.rmin) + P.rmin;
        R(lpc,:) = x;
        phaseVec = x(1)*P.dataX + x(2)*P.dataXSq + x(3)*P.dataXCb;
        qc = sin(2*pi*phaseVec);
        sampFreq = 1 / (P.dataX(2) - P.dataX(1));
        [qc_normalized, ~] = normsig4psd(qc, sampFreq, P.psdPosFreq, 1);
        llr = innerprodpsd(P.dataY, qc_normalized, sampFreq, P.psdPosFreq);
        F(lpc) = -llr^2;
    end
end

function innProd = innerprodpsd(xVec, yVec, sampFreq, psdVals)
    nSamples = length(xVec);
    kNyq = floor(nSamples/2)+1;
    fftX = fft(xVec);
    fftY = fft(yVec);
    negFStrt = 1-mod(nSamples,2);
    psdVec4Norm = [psdVals, psdVals((kNyq-negFStrt):-1:2)];
    dataLen = sampFreq*nSamples;
    innProd = (1/dataLen)*(fftX./psdVec4Norm)*fftY';
    innProd = real(innProd);
end
