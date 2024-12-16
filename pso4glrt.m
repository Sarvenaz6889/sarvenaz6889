function outStruct = pso4glrt(inParams, options, numRuns)
% PSO4GLRT: Particle Swarm Optimization for GLRT Quadratic Chirp Estimation
% 
% INPUTS:
%   inParams - Struct with input parameters:
%       .dataX      - Time vector
%       .dataY      - Data vector
%       .rmin       - Lower bounds for parameters [a1, a2, a3]
%       .rmax       - Upper bounds for parameters [a1, a2, a3]
%   options  - Struct with optimization settings:
%       .maxSteps   - Maximum number of PSO iterations
%   numRuns  - Number of PSO runs
%
% OUTPUTS:
%   outStruct - Struct with optimization results:
%       .bestQcCoefs   - Best estimated coefficients [a1, a2, a3]
%       .bestFitness   - Best fitness value
%       .bestSig       - Best estimated signal
%       .allRunsOutput - All PSO run results (positions, fitness values)

    % --- Extract Inputs ---
    timeVec = inParams.dataX;    % Time vector
    dataVec = inParams.dataY;    % Data vector
    rmin = inParams.rmin;        % Lower bounds
    rmax = inParams.rmax;        % Upper bounds
    maxSteps = options.maxSteps; % Max PSO steps
    
    % --- PSO Parameters ---
    numParticles = 30;           % Number of particles in the swarm
    numParams = 3;               % Number of parameters to optimize (a1, a2, a3)
    
    % Initialize output structure
    allRunsOutput = struct([]);
    
    % --- Run PSO Optimization Multiple Times ---
    for i = 1:numRuns
        disp(['Starting PSO Run ', num2str(i), '...']);
        
        % Call the PSO function
        [bestPos, bestFitness] = simplePSO(@glrtFitness, numParams, rmin, rmax, maxSteps, numParticles, timeVec, dataVec);
        
        % Store results for this run
        allRunsOutput(i).bestPos = bestPos;
        allRunsOutput(i).bestFitness = bestFitness;
        disp(['Run ', num2str(i), ' Best Fitness: ', num2str(bestFitness)]);
    end
    
    % --- Identify the Best Run ---
    fitnessValues = [allRunsOutput.bestFitness];
    [~, bestRunIdx] = min(fitnessValues); % Find the index of the best run
    
    % Best parameters and signal
    bestQcCoefs = allRunsOutput(bestRunIdx).bestPos;
    bestSig = glrtSignal(bestQcCoefs, timeVec);
    
    % --- Populate Output Struct ---
    outStruct.bestQcCoefs = bestQcCoefs;       % Best parameters
    outStruct.bestFitness = fitnessValues(bestRunIdx); % Best fitness value
    outStruct.bestSig = bestSig;               % Best estimated signal
    outStruct.allRunsOutput = allRunsOutput;   % All runs' results
end
