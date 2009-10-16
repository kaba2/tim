% TRANSFER_ENTROPY_T
% A temporal transfer entropy estimate from samples.
%
% I = transfer_entropy_t(X, Y, W, 
%       timeWindowRadius, xLag, yLag, wLag, k, threads)
%
% where
%
% X, Y, and W are cell-arrays of arbitrary dimension whose linearizations
% contain q trials of the signals X, Y, and W, respectively. 
%
% TIMEWINDOWRADIUS determines the radius of the time-window in samples 
% inside which samples are taken into consideration to the estimate at 
% time instant t. This allows the estimate to be adaptive to temporal changes.
% If no such changes should happen, better accuracy can be 
% achieved by either setting 'timeWindowRadius' maximally wide
% or by using the transfer_entropy() function instead.
%
% XLAG, YLAG, and WLAG are the lags in samples applied to 
% signals X, Y, and W, respectively.
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
% m-dimensional signal. The signals contained in X (or Y or W) 
% must all have equal dimensionality, but their number of samples may vary. 
% If the number of samples varies with trials, the function uses 
% the minimum sample count among the trials of X, Y, and W.
% The number of trials in X, Y, and W must be equal.

% Description: Temporal transfer entropy estimation
% Documentation: tim_matlab.txt

function I = transfer_entropy_t(X, Y, W, ...
    timeWindowRadius, xLag, yLag, wLag, k, threads)

if nargin < 4
    error('Not enough input arguments.');
end

if ~iscell(X) || ~iscell(Y) || ~iscell(W)
    error('X, Y, or W is not a cell-array.');
end

if nargin >= 5 && nargin < 7
    error('Lags must be specified all at once to avoid errors.');
end

if nargin < 5
    xLag = 0;
    yLag = 0;
    wLag = 0;
end

if nargin < 8
    k = 1;
end

if nargin < 9
    threads = 1;
end

I = entropy_combination_t(...
    [W(:), X(:), Y(:)]', ...
    [1, 2, 1; 2, 3, 1; 2, 2, -1], ...
    timeWindowRadius, ...
    [wLag, xLag, yLag], k, threads);
