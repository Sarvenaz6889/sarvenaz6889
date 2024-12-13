% Set up input parameters
nSamples = 512;
sampFreq = 512; % Hz
SNR = 10; 
a1 = 10; 
a2 = 3; 
a3 = 3;

% Search range of phase coefficients
rmin = [1, 1, 1];
rmax = [180, 10, 10];
dataX = (0:(nSamples-1))/Fs;
dataLen = nSamples/Fs;
sigVec = qcsigfunc(dataX,SNR,[a1,a2,a3]);
[sigVec, ~] = normsig4psd(sigVec, Fs, psdPosFreq, 10);
kNyq = floor(nSamples/2)+1;
posFreq = (0:(kNyq-1))*(1/dataLen);
noisePSD = @(f) (f>=50 & f<=100).*(f-50).*(100-f)/625 + 1;
psdPosFreq = noisePSD(posFreq);
dataY = sigVec+outnoise;
[dataY, sig] = crcbgenqcdata(dataX,snr,[a1,a2,a3]);

fltrOrder = 500;
outnoise = statgaussnoisegen(nSamples,[posFreq(:),psdPosFreq(:)],fltrOrder,Fs);

sigVec = qcsigfunc(dataX,SNR,[a1,a2,a3]);
[sigVec, ~] = normsig4psd(sigVec, Fs, psdPosFreq, 10);

% Number of independent PSO runs
nRuns = 8;

rng('default');

% Input parameters 
inParams = struct('dataX', dataX,...
                  'dataY', dataY,...
                  'dataXSq',dataX.^2,...
                  'dataXCb',dataX.^3,...
                  'rmin',rmin,...
                  'rmax',rmax);
% CRCBQCHRPPSO runs PSO on the CRCBQCHRPFITFUNC fitness function. As an
% illustration of usage, we change one of the PSO parameters from its
% default value.
outStruct = crcbqcpso(inParams,struct('maxSteps',2000),nRuns);

% Plots
figure;
hold on;
plot(dataX,dataY,'.');
plot(dataX,sig);
for lpruns = 1:nRuns
    plot(dataX,outStruct.allRunsOutput(lpruns).estSig,'Color',[51,255,153]/255,'LineWidth',4.0);
end
plot(dataX,outStruct.bestSig,'Color',[76,153,0]/255,'LineWidth',2.0);
legend('Data','Signal',...
       ['Estimated signal: ',num2str(nRuns),' runs'],...
       'Estimated signal: Best run');
disp(['Estimated parameters: a1=',num2str(outStruct.bestQcCoefs(1)),...
                             '; a2=',num2str(outStruct.bestQcCoefs(2)),...
                             '; a3=',num2str(outStruct.bestQcCoefs(3))]);
