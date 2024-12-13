% Define the problem's dimensions and bounds
nDim = 30;
rmin = -100;
rmax = 100;

% Fitness function parameters
fitFuncParams = struct('rmin', -100 * ones(1, nDim), ...
                       'rmax', 100 * ones(1, nDim));

% Fitness function handle
fitFuncHandle = @(x) crcbpsotestfunc_sphere(x, fitFuncParams);

% PSO parameters
psoParams = struct('popSize', 40, 'maxSteps', 2000, 'c1', 2, 'c2', 2, ...
                   'maxVelocity', 0.5, 'startInertia', 0.9, 'endInertia', 0.4, ...
                   'boundaryCond', '', 'nbrhdSz', 3);

% Output level (verbosity)
outputLvl = 2;

% Run PSO to minimize the Sphere function
psoOutput = crcbpso(fitFuncHandle, nDim, psoParams, outputLvl);

% Display results
disp('Global Minimum Found:');
disp(psoOutput.bestFitness);

disp('Minimizer Coordinates (Standardized):');
disp(psoOutput.bestLocation);

% Generate a grid for plotting
x = linspace(-10, 10, 100); % Search range in X-axis
y = linspace(-10, 10, 100); % Search range in Y-axis
[X, Y] = meshgrid(x, y);

% Standardize coordinates for plotting
Xstd = (X - fitFuncParams.rmin(1)) / (fitFuncParams.rmax(1) - fitFuncParams.rmin(1));
Ystd = (Y - fitFuncParams.rmin(1)) / (fitFuncParams.rmax(1) - fitFuncParams.rmin(1));

% Flatten Xstd and Ystd into column vectors
coords = [Xstd(:), Ystd(:)];

% Evaluate the fitness function at each grid point
Z = arrayfun(@(i) fitFuncHandle(coords(i, :)), 1:size(coords, 1));

% Reshape the fitness values back into the grid shape for plotting
Z = reshape(Z, size(X));

% Create surface plot
figure;
surf(X, Y, Z);
title('3D Surface Plot of Fitness Function');
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Fitness Value');
shading interp; % Smooth the surface
