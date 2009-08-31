% TEMPORAL_MUTUAL_INFORMATION 
% A temporal mutual information estimate from samples.
%
% I = temporal_mutual_information(X, Y, timeWindowRadius, yLag, k, threads)
%
% where
%
% X is an arbitrary-dimensional cell-array whose linearization
% contains q trials of signal x.
%
% Y is an arbitrary-dimensional cell-array whose linearization 
% contains q trials of signal y.
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
% YLAG is the lag in samples which is applied to signal Y.
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
% m-dimensional signal. The dimension of X and Y need not coincide.
% However, the number of trials has to coincide.
% If the number of samples varies with trials, the function uses 
% the minimum sample count among the trials of X and Y.

% Description: Temporal mutual information estimation
% Documentation: tim_matlab.txt

function I = temporal_mutual_information(X, Y, timeWindowRadius, ...
    yLag, k, threads)

% The limit for the dimension is arbitrary, but
% protects for the case when the user accidentally
% passes the transpose of the intended data.
maxDimension = 32;

if nargin < 3
    error('Not enough input arguments.');
end

if nargin > 6
    error('Too many input arguments.');
end

if nargout > 1
    error('Too many output arguments.');
end

if nargin < 4
    yLag = 0;
end

if nargin < 5
    k = 1;
end

if nargin < 6
    threads = 1;
end

if ~iscell(X)
    error('X must be a cell array.');
end

if ~iscell(Y)
    error('Y must be a cell array.');
end

signals = prod(size(X));

if prod(size(Y)) ~= signals
    error('X and Y must contain the same number of signals');
end

for i = 1 : signals
    if ~isa(X{i}, 'double')
        error('Some signal of X is not of type double.');
    end

    if size(X{i}, 1) > maxDimension
        error(['Some signal of X has dimension greater than ', ...
            int2str(maxDimension), '.']);
    end
end

for i = 1 : signals
    if ~isa(Y{i}, 'double')
        error('Some signal of Y is not of type double.');
    end

    if size(Y{i}, 1) > maxDimension
        error(['Some signal of Y has dimension greater than ', ...
            int2str(maxDimension), '.']);
    end
end

if size(timeWindowRadius, 1) ~= 1 || ...
   size(timeWindowRadius, 2) ~= 1
    error('TIMEWINDOWRADIUS must be a scalar integer.');
end

if timeWindowRadius < 0
    error('TIMEWINDOWRADIUS must be non-negative.');
end

if size(yLag, 1) ~= 1 || ...
   size(yLag, 2) ~= 1
    error('YLAG must be a scalar integer.');
end

if yLag < 0
    error('YLAG must be non-negative.');
end

if size(k, 1) ~= 1 || ...
   size(k, 2) ~= 1
    error('K must be a scalar integer.');
end

if k < 1
    error('K must be at least 1.');
end

if size(threads, 1) ~= 1 || ...
   size(threads, 2) ~= 1
    error('THREADS must be a scalar integer.');
end

if threads < 1
    error('THREADS must be at least 1.');
end

I = timTemporalMutualInformation(X, Y, timeWindowRadius, ...
    yLag, k, threads);

