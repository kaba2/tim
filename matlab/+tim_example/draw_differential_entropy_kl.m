% DRAW_DIFFERENTIAL_ENTROPY_KL_DISTRIBUTION
% Distribution of differential_entropy_kl -estimates.
%
% draw_differential_entropy_kl('key', value, ...)
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
% D ('d') is a positive integer which denotes the dimension of the point-set.
% Default: 4
%
% COV ('cov') is a (DxD) positive-definite real matrix which specifies the
% covariance matrix of the point-set.
% Default: eye(d, d)
%
% M ('m') is a non-negative integer which denotes the number of estimates 
% to compute. The histogram is formed from these estimates.
% Default: 1000
%
% As a trial, we generate a point-set from a multi-variate normal 
% distribution with identity covariance, and compute a differential 
% entropy estimate via differential_entropy_kl. We will then accumulate 
% many such trials to numerically reveal the distribution of the 
% estimator for the chosen arguments.

% Description: Draws the distribution of differential entropy estimates

function draw_differential_entropy_kl(varargin)

% Optional input arguments
k = 1;
d = 4;
n = 1000;
m = 1000;
cov = {};
eval(tim.process_options({'k', 'd', 'n', 'm', 'cov'}, varargin));

if iscell(cov)
    cov = eye(d, d);
end

pastelsys.concept_check(...
    k, 'integer', ...
    k, 'positive', ...
    d, 'integer', ...
    d, 'positive', ...
    n, 'integer', ...
    n, 'non_negative', ...
    m, 'integer', ...
    m, 'non_negative');

% Accumulate the estimates for each trial.
hSet = zeros(1, m);
for i = 1 : m
    % Generate a point-set randomly from the multi-variate
    % normal distribution.
	A = tim.random_normal(d, n, 'cov', cov);
    % Compute the differential entropy estimate.
	hSet(i) = tim.differential_entropy_kl(A, 'k', k);
end

% Compute the sample mean of the estimates.
hMean = mean(hSet);

% Compute the standard deviation of the estimates.
hDeviation = std(hSet);

% Compute the correct differential entropy; the
% case for the multi-variate normal distribution can be
% solved analytically.
hCorrect = tim.differential_entropy_normal(d, 'detCov', det(cov));

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
title(['Distribution of the KL differential entropy -estimate', ...
    ', k = ', num2str(k), ...
    ', d = ', num2str(d), ...
    ', n = ', num2str(n)]);
xlabel(['KL differential entropy -estimate', ...
    ', mean = ', num2str(hMean), ...
    ', deviation = ', num2str(hDeviation)]);
ylabel('Samples');
legend('Samples', 'Sample mean', 'Correct value');
axis([xMin, xMax, yLimits(1), yLimits(2)]);
hold off;