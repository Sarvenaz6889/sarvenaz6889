function glrtValue = glrtqcsig(dataVec, timeVec, psdPosFreq, qcCoefs)
    % Generate the unit norm signal (template)
    sigVec = crcbgenqcsig(timeVec, 1, qcCoefs);
    
    % Calculate sampling frequency
    sampFreq = 1 / (timeVec(2) - timeVec(1));
    
    % Normalize the template vector to unit norm
    [templateVec, ~] = normsig4psd(sigVec, sampFreq, psdPosFreq, 1);
    
    % Calculate inner product of the data with the unit norm template
    llr = innerprodpsd(dataVec, templateVec, sampFreq, psdPosFreq);
    
    % GLRT is the square of the inner product
    glrtValue = llr^2;
end