% MUTUAL_INFORMATION_PT
% A temporal partial mutual information estimate from samples.
%
% I = mutual_information_pt(
%         X, Y, Z, timeWindowRadius, 
%         xLag, yLag, zLag, k, threads)
%
% where
%
% X, Y, and Z are arbitrary-dimensional cell-arrays whose 
% linearizations contain q trials of signal x, y, and z, 
% respectively.
%
% TIMEWINDOWRADIUS determines the radius of the time-window in samples 
% inside which samples are taken into consideration to the estimate at 
% time instant t. This allows the estimate to be adaptive to temporal changes.
% If no such changes should happen, better accuracy can be 
% achieved by either setting 'timeWindowRadius' maximally wide
% or by using the mutual_information_p() function instead.
%
% XLAG, YLAG and ZLAG are the lags in samples applied to signals
% X, Y and Z, respectively.
%
% K determines which k:th nearest neighbor the algorithm
% uses for estimation. Default 1.
%
% THREADS determines the number of threads to use for parallelization.
% To fully take advantage of multiple cores in your machine, set this
% to the number of cores in your machine. Note however that this makes 
% your computer unresponsive to other tasks. When you need responsiveness, 
% spare one core for other work. Default maxNumCompThreads.
%
% Each signal is a real (m x n)-matrix that contains n samples of an
% m-dimensional signal. The signals contained in X (or Y or Z) must all 
% have equal dimensionality, but their number of samples may vary. 
% If the number of samples varies with trials, the function uses 
% the minimum sample count among the trials of X, Y, and Z.
% The number of trials in X, Y, and Z must be equal.

% Description: Temporal partial mutual information estimation
% Documentation: tim_matlab_matlab.txt

function I = mutual_information_pt(...
    X, Y, Z, timeWindowRadius, xLag, yLag, zLag, k, threads)

if nargin < 4
    error('Not enough input arguments.');
end

if nargin >= 5 && nargin < 7
    error('Lags must be specified all at once to avoid errors.');
end

if nargin < 5
    xLag = 0;
    yLag = 0;
    zLag = 0;
end

if nargin < 8
    k = 1;
end

if nargin < 9
    threads = maxNumCompThreads;
end

if ~iscell(X) || ~iscell(Y) || ~iscell(Z)
    error('X, Y, or Z is not a cell-array.');
end

if numel(X) ~= numel(Y) || numel(X) ~= numel(Z)
    error('The number of trials in X, Y, and Z differ.');
end

% Pass parameter error checking to entropy_combination.

I = entropy_combination_t(...
    [X(:), Z(:), Y(:)]', ...
    [1, 2, 1; 2, 3, 1; 2, 2, -1], timeWindowRadius, ...
    [xLag, zLag, yLag], k, threads);
