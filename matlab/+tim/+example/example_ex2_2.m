% Non-linear time series analysis
% Exercise 2.2

% Diverging from the exercise, I'll map to the [-1, 1] range instead 
% of [0, 1], so that the constant-component will not distract in the 
% power spectrum.

clear all;
close all;

% Uniform random
% --------------

n = 4096;
dataSet = 2 * rand(1, n) - 1;
h = tim.example.visualize(dataSet);
fprintf('Uniform random\n');
fprintf('Sample mean: %0.4f\n', mean(dataSet));
fprintf('Sample standard deviation: %0.4f\n', std(dataSet));
fprintf('Mean: %0.4f\n', 0);
fprintf('Standard deviation: %0.4f\n\n', sqrt(1 / 3));

% Distorted Ulam map
% ------------------

dataSet = zeros(1, n);
dataSet(1) = 0.1;
for i = 2 : n
    % This the Ulam map.
    dataSet(i) = 1 - 2 * dataSet(i - 1)^2;
end

% The data is measured through a non-linear
% measurement function.
for i = 1 : n
    dataSet(i) = 2 * (acos(-dataSet(i)) / pi) - 1;
end

tim.example.visualize(dataSet);
fprintf('Distorted Ulam map\n');
fprintf('Sample mean: %0.4f\n', mean(dataSet));
fprintf('Sample standard deviation: %0.4f\n', std(dataSet));

% So the distorted Ulam map looks very much like uniform
% random noise.