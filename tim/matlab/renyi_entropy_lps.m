% RENYI_ENTROPY_LPS
% A Renyi entropy estimate from samples
% using Leonenko-Pronzato-Savani nearest neighbor estimator.
%
% H = renyi_entropy_lps(S, q, epsilon, k, threads)
%
% where
%
% S is an arbitrary-dimensional cell-array whose linearization contains
% q trials of a signal. Each signal is a real (m x n)-matrix that 
% contains n samples of an m-dimensional signal.
%
% Q is the power in the definition Renyi entropy.
% In case Q = 1, differential_entropy_kl() is used to
% compute the result instead. Default 2.
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
% spare one core for other work. Default maxNumCompThreads.

% Description: Renyi entropy estimation
% Detail: Leonenko-Pronzato-Savani nearest neighbor estimator
% Documentation: tim_matlab.txt

function H = renyi_entropy_lps(S, q, epsilon, k, threads)

if nargin < 1
    error('Not enough input arguments.');
end

if nargin > 5
    error('Too many input arguments.');
end

if nargout > 1
    error('Too many output arguments.');
end

if nargin < 2
	q = 2
end

if nargin < 3
    epsilon = 0;
end

if nargin < 4
    k = 1;
end

if nargin < 5
    threads = maxNumCompThreads;
end

checkSignalSet(S);

if size(q, 1) ~= 1 || ...
   size(q, 2) ~= 1
    error('Q must be a scalar.');
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
    error('EPSILON must be non-negative.');
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

H = timRenyiEntropyLps(S, q, epsilon, k, threads);
