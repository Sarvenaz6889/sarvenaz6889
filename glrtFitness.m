function fitness = glrtFitness(params, timeVec, dataVec)
    % Calculate the signal using params
    estSig = glrtSignal(params, timeVec);
    % Mean squared error between data and estimated signal
    fitness = sum((dataVec - estSig).^2);
end
