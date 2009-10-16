% TRANSFER_ENTROPY_PT
% A temporal partial transfer entropy estimate from samples.
%
% I = transfer_entropy_pt(X, Y, Z, W, 
%       timeWindowRadius, xLag, yLag, zLag, wLag, k, threads)
%
% where
%
% X, Y, Z, and W are cell-arrays of arbitrary dimension whose linearization
% contains q trials of the signals X, Y, Z, and W, respectively. 
%
% TIMEWINDOWRADIUS determines the radius of the time-window in samples 
% inside which samples are taken into consideration to the estimate at 
% time instant t. This allows the estimate to be adaptive to temporal 
% changes.
% If no such changes should happen, better accuracy can be 
% achieved by either setting 'timeWindowRadius' maximally wide
% or by using the transfer_entropy() function instead.
%
% XLAG, YLAG, ZLAG, and WLAG are the lags in samples applied to 
% signals X, Y, Z, and W, respectively.
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
% m-dimensional signal. The signals contained in X (or Y or Z or W) 
% must all have equal dimensionality, but their number of samples may vary. 
% If the number of samples varies with trials, the function uses 
% the minimum sample count among the trials of X, Y, Z, and W.
% The number of trials in X, Y, Z, and W must be equal.

% Description: Temporal partial transfer entropy estimation
% Documentation: tim_matlab.txt

function I = transfer_entropy_pt(X, Y, Z, W, ...
    timeWindowRadius, xLag, yLag, zLag, wLag, k, threads)

if nargin < 5
    error('Not enough input arguments.');
end

if ~iscell(X) || ~iscell(Y) || ~iscell(Z) || ~iscell(W)
    error('X, Y, Z, or W is not a cell-array.');
end

if nargin >= 6 && nargin < 9
    error('Lags must be specified all at once to avoid errors.');
end

if nargin < 6
    xLag = 0;
    yLag = 0;
    zLag = 0;
    wLag = 0;
end

if nargin < 10
    k = 1;
end

if nargin < 11
    threads = 1;
end

% Pass parameter error checking to entropy_combination.

I = entropy_combination_t(...
    [W(:), X(:), Z(:), Y(:)]', ...
    [1, 3, 1; 2, 4, 1; 2, 3, -1], ...
    timeWindowRadius, ...
    [wLag, xLag, zLag, yLag], k, threads);
