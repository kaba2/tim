% DIFFERENTIAL_ENTROPY_KL
% A differential entropy estimate from samples
% using Kozachenko-Leonenko nearest neighbor estimator.
%
% H = differential_entropy_kl(S, epsilon, k, threads)
%
% where
%
% S is an arbitrary dimensional cell array whose linearization contains
% q trials of a signal. Each signal is a real (m x n)-matrix that 
% contains n samples of an m-dimensional signal.
%
% EPSILON is the maximum relative error in distance that
% nearest neighbor searching is allowed to result in.
% Higher tolerances result in enhanced performance, but
% increases errors in the estimate. Default 0.
%
% K determines which k:th nearest neighbor the algorithm
% uses for estimation. Default 1.
%
% THREADS determines the number of threads to use for parallelization.
% To fully take advantage of multiple cores in your machine, set this
% to the number of cores in your machine. Note however that this makes 
% your computer unresponsive to other tasks. When you need responsiveness, 
% spare one core for other work. Default 1 (no parallelization).

% Description: Differential entropy estimation
% Detail: Kozachenko-Leonenko nearest neighbor estimator
% Documentation: tim_matlab.txt

function H = differential_entropy_kl(S, epsilon, k, threads)

% The limit for the dimension is arbitrary, but
% protects for the case when the user accidentally
% passes the transpose of the intended data.
maxDimension = 32;

if nargin < 1
    error('Not enough input arguments.');
end

if nargin > 4
    error('Too many input arguments.');
end

if nargout > 1
    error('Too many output arguments.');
end

if nargin < 2
    epsilon = 0;
end

if nargin < 3
    k = 1;
end

if nargin < 4
    threads = 1;
end

if ~iscell(S)
    error('S must be a cell array.');
end

signals = prod(size(S));

if signals < 1
    error('S must contain at least 1 signal.');
end

dimension = size(S{1}, 1);

for i = 1 : signals
    if ~isa(S{i}, 'double')
        error('Some signal of S is not of type double.');
    end

    if size(S{i}, 1) ~= dimension
        error(['The dimensions of the trials do not match.']);
    end
    
    if size(S{i}, 1) > maxDimension
        error(['Some signal of S has dimension greater than ', ...
            int2str(maxDimension), '.']);
    end
end

if size(epsilon, 1) ~= 1 || ...
   size(epsilon, 2) ~= 1
    error('EPSILON must be a scalar.');
end

if size(k, 1) ~= 1 || ...
   size(k, 2) ~= 1
    error('K must be a scalar integer.');
end

if epsilon < 0
    error('EPSILON must be non-negative'.');
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

H = timDifferentialEntropyKl(S, epsilon, k, threads);
