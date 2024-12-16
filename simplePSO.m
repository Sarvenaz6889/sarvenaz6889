function [bestPos, bestFitness] = simplePSO(fitnessFunc, numParams, rmin, rmax, maxSteps, numParticles, timeVec, dataVec)
    % SIMPLEPSO: A basic Particle Swarm Optimization implementation
    % Initialize positions and velocities
    positions = rand(numParticles, numParams) .* (rmax - rmin) + rmin;
    velocities = zeros(size(positions));
    personalBestPos = positions;
    personalBestFitness = inf(numParticles, 1);
    [globalBestFitness, bestIdx] = min(personalBestFitness);
    globalBestPos = personalBestPos(bestIdx, :);

    % PSO loop
    for step = 1:maxSteps
        for i = 1:numParticles
            % Evaluate fitness
            fitness = fitnessFunc(positions(i, :), timeVec, dataVec);
            if fitness < personalBestFitness(i)
                personalBestFitness(i) = fitness;
                personalBestPos(i, :) = positions(i, :);
            end
            if fitness < globalBestFitness
                globalBestFitness = fitness;
                globalBestPos = positions(i, :);
            end
        end
        % Update velocities and positions
        w = 0.7; c1 = 1.5; c2 = 1.5; % PSO parameters
        r1 = rand(numParticles, numParams);
        r2 = rand(numParticles, numParams);
        velocities = w * velocities ...
                     + c1 * r1 .* (personalBestPos - positions) ...
                     + c2 * r2 .* (globalBestPos - positions);
        positions = positions + velocities;
    end

    % Return the best position and fitness
    bestPos = globalBestPos;
    bestFitness = globalBestFitness;
end
