% MUTUAL_INFORMATION_T
% A temporal mutual information estimate from samples.
%
% I = mutual_information_t(X, Y, timeWindowRadius, 
%       xLag, yLag, k, threads)
%
% where
%
% X and Y are arbitrary-dimensional cell-arrays whose linearizations
% contain q trials of signal X and Y, respectively.
%
% TIMEWINDOWRADIUS determines the radius of the time-window inside which
% samples are taken into consideration to the mutual information
% estimate at time instant t. The time window at time instant t
% is given by [t - timeWindowRadius, t + timeWindowRadius]. 
% This allows the estimate to be adaptive to temporal changes in mutual 
% information. If no such changes should happen, better accuracy can be 
% achieved by either setting 'timeWindowRadius' to the number of samples 
% or using the mutual_information() function instead.
%
% XLAG and YLAG are the lags in samples which are applied 
% to signals X and Y.
%
% K determines which k:th nearest neighbor the algorithm
% uses for estimation. Default 1.
%
% THREADS determines the number of threads to use for parallelization.
% To fully take advantage of multiple cores in your machine, set this
% to the number of cores in your machine. Note however that this makes 
% your computer unresponsive to other tasks. When you need responsiveness, 
% spare one core for other work. Default 1 (no parallelization).
%
% Each signal is a real (m x n)-matrix that contains n samples of an
% m-dimensional signal. The signals contained in X (Y) must all have equal
% dimensionality, but their number of samples may vary. 
% If the number of samples varies with trials, the function uses 
% the minimum sample count among the trials of X and Y.
% The number of trials in X and Y must be equal.

% Description: Temporal mutual information estimation
% Documentation: tim_matlab.txt

function I = mutual_information_t(X, Y, timeWindowRadius, ...
    xLag, yLag, k, threads)

if nargin < 3
    error('Not enough input arguments.');
end

if ~iscell(X) || ~iscell(Y)
    error('X or Y is not a cell-array.');
end

if nargin >= 4 && nargin < 5
    error('Lags must be specified all at once to avoid errors.');
end

if nargin < 4
    xLag = 0;
    yLag = 0;
end

if nargin < 6
    k = 1;
end

if nargin < 7
    threads = 1;
end

I = entropy_combination_t(...
    [X(:), Y(:)]', ...
    [1, 1, 1; 2, 2, 1], timeWindowRadius, ...
    [xLag, yLag], k, threads);
