% DRAW_MUTUAL_INFORMATION
% Draws the distribution of mutual_information -estimates.
%
% draw_mutual_information('key', value, ...)
%
% Optional input arguments in 'key'-value pairs:
%
% K ('k') is a positive integer which denotes the number of nearest 
% neighbors to be used by the estimator. Greater number decreases bias,
% but increases variance.
% Default: 1
%
% N ('n') is a non-negative integer which denotes the number 
% of points to generate.
% Default: 1000
%
% M ('m') is a non-negative integer which denotes the number of estimates 
% to compute. The histogram is formed from these estimates.
% Default: 1000
%
% COV ('cov') is a positive-definite (D x D) real matrix which
% denotes the covariance matrix of the random variable (X, Y).
% The D also determines the dimension of the random variable (X, Y).
% Default: eye(4, 4)
%
% DX ('dx') is a positive integer which denotes the dimension of 
% the random variable X.
% Default: floor(D / 2)
%
% As a trial, we generate point-sets X and Y from a multi-variate normal 
% distribution with identity covariance, and compute a differential 
% entropy estimate via differential_entropy_kl. We will then accumulate 
% many such trials to numerically reveal the distribution of the 
% estimator for the chosen arguments.

% Description: Draws the distribution of mutual_information -estimates

function draw_mutual_information(varargin)

% Optional input arguments
k = 1;
n = 1000;
m = 1000;
cov = eye(4, 4);
dx = {};
eval(tim.process_options({'k', 'n', 'm', 'cov', 'dx'}, varargin));

d = size(cov, 1);
if iscell(dx)
    dx = floor(d / 2);
end

pastelmatlab.concept_check(...
    k, 'integer', ...
    k, 'positive', ...
    n, 'integer', ...
    n, 'non_negative', ...
    m, 'integer', ...
    m, 'non_negative', ...
    dx, 'integer', ...
    dx, 'positive', ...
    cov, 'real_matrix', ...
    cov, 'square_matrix');

% Accumulate the estimates for each trial.
hSet = zeros(1, m);
for i = 1 : m
    % Generate a point-set randomly from the multi-variate
    % normal distribution.
	XY = tim.random_normal(d, n, 'cov', cov);
    X = XY(1 : dx, :);
    Y = XY((dx + 1) : end, :);
    % Compute the mutual information estimate.
	hSet(i) = tim.mutual_information(X, Y, 'k', k);
end

% Compute the sample mean of the estimates.
hMean = mean(hSet);

% Compute the standard deviation of the estimates.
hDeviation = std(hSet);

xCov = cov(1 : dx, 1 : dx);
yCov = cov((dx + 1) : end, (dx + 1) : end);

% Compute the correct mutual information; the
% case for the multi-variate normal distribution can be
% solved analytically.
hCorrect = tim.mutual_information_normal(det(xCov) * det(yCov), det(cov));

xMin = hCorrect - 0.25;
xMax = hCorrect + 0.25;

% Draw a histogram of the estimates.
hist(hSet, 100);
hold on
yLimits = [0, 50];

% Draw the sample mean as a vertical line.
line([hMean, hMean], yLimits, 'Color', [1, 0, 0]);
% Draw the correct value as a vertical line.
line([hCorrect, hCorrect], yLimits, 'Color', [0, 1, 0]);

% Name the axes and the figure.
title(['Distribution of the mutual information -estimate', ...
    ', k = ', num2str(k), ...
    ', d = ', num2str(d), ...
    ', n = ', num2str(n)]);
xlabel(['mutual information -estimate', ...
    ', mean = ', num2str(hMean), ...
    ', deviation = ', num2str(hDeviation)]);
ylabel('Samples');
legend('Samples', 'Sample mean', 'Correct value');
axis([xMin, xMax, yLimits(1), yLimits(2)]);
hold off;