% TEMPORAL_PARTIAL_MUTUAL_INFORMATION 
% A temporal partial mutual information estimate from samples.
%
% I = temporal_partial_mutual_information(
%         X, Y, Z, timeWindowRadius, yLag, zLag, k, threads)
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
% or by using the partial_mutual_information() function instead.
%
% YLAG and ZLAG are the lags in samples applied to signals
% y and z, respectively.
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
% m-dimensional signal. The dimensions of X, Y, and Z need not coincide.
% However, the number of trials has to coincide.
% If the number of samples varies with trials, the function uses 
% the minimum sample count among the trials of X and Y.

function I = temporal_partial_mutual_information(...
    X, Y, Z, timeWindowRadius, yLag, zLag, k, threads)

% The limit for the dimension is arbitrary, but
% protects for the case when the user accidentally
% passes the transpose of the intended data.
maxDimension = 32;

if nargin < 4
    error('Not enough input arguments.');
end

if nargin > 8
    error('Too many input arguments.');
end

if nargout > 1
    error('Too many output arguments.');
end

if nargin < 5
    yLag = 0;
end

if nargin < 6
    zLag = 0;
end

if nargin < 7
    k = 1;
end

if nargin < 8
    threads = 1;
end

if ~iscell(X)
    error('X must be a cell array.');
end

if ~iscell(Y)
    error('Y must be a cell array.');
end

if ~iscell(Z)
    error('Z must be a cell array.');
end

signals = prod(size(X));

if prod(size(Y)) ~= signals || prod(size(Z)) ~= signals
    error('X, Y, and Z must contain the same number of signals.');
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

for i = 1 : signals
    if ~isa(Z{i}, 'double')
        error('Some signal of Z is not of type double.');
    end

    if size(Z{i}, 1) > maxDimension
        error(['Some signal of Z has dimension greater than ', ...
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

if size(zLag, 1) ~= 1 || ...
   size(zLag, 2) ~= 1
    error('ZLAG must be a scalar integer.');
end

if zLag < 0
    error('ZLAG must be non-negative.');
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

%I = timTemporalPartialMutualInformation(...
%    X, Y, Z, timeWindowRadius, yLag, zLag, k, threads);

