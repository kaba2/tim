% MUTUAL_INFORMATION A mutual information estimate from samples.
% I = mutual_information(S, k, threads)
% where
% S is a cell array of arbitrary dimension that contains p signals 
% (it is addressed as a 1d cell-array).
% Each signal is a real (m x n)-matrix that contains n samples of an
% m-dimensional signal. If the number of samples varies with each
% signal, the function uses the minimum sample count among the signals.
% K determines which k:th nearest neighbor the algorithm
% uses for estimation. Default 1.
% THREADS determines the number of threads to use for parallelization.
% To fully take advantage of multiple cores in your machine, set this
% to the number of cores in your machine. Note however that this makes 
% your computer unresponsive to other tasks. When you need responsiveness, 
% spare one core for other work. Default 1 (no parallelization).

function I = mutual_information(S, k, threads)

if nargin < 1
    error('Not enough input arguments.');
end

if nargin > 3
    error('Too many input arguments.');
end

if nargout > 1
    error('Too many output arguments.');
end

if nargin < 2
    k = 1;
end

if nargin < 3
    threads = 1;
end

if ~iscell(S)
    error('S must be a cell array.');
end

signals = prod(size(S));

if signals < 2
    error('S must contain at least 2 signals.');
end

for i = 1 : signals
    if ~isa(S{i}, 'double')
        error('A signal must be of type double.');
    end

    % The limit for the dimension is arbitrary, but
    % protects for the case when the user accidentally
    % passes the transpose of the intended data.
    if size(S{i}, 1) > 32
        error('A signal must have dimension less than 32.');
    end
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

I = timMutualInformation(S, k, threads);
