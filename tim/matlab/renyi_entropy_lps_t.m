% RENYI_ENTROPY_LPS_T
% A temporal Renyi entropy estimate from samples
% using Leonenko-Pronzato-Savani nearest neighbor estimator.
%
% H = renyi_entropy_kl_t(
%     S, timeWindowRadius, q, epsilon, kSuggestion, threads)
%
% where
%
% S is an arbitrary-dimensional cell-array whose linearization contains
% q trials of a signal. Each signal is a real (m x n)-matrix that 
% contains n samples of an m-dimensional signal.
%
% Q is the power in the definition Renyi entropy.
% If Q = 1, differential_entropy_kl_t() is used to
% compute the result instead. 
% If Q < 1, there are huge errors in the estimation.
% Default 2.
%
% TIMEWINDOWRADIUS determines the radius of the time-window in samples 
% inside which samples are taken into consideration to the estimate at 
% time instant t. This allows the estimate to be adaptive to temporal changes.
% If no such changes should happen, better accuracy can be 
% achieved by either setting 'timeWindowRadius' maximally wide
% or by using the renyi_entropy() function instead.
%
% EPSILON is the maximum relative error in distance that
% nearest neighbor searching is allowed to result in.
% Higher tolerances result in enhanced performance, but
% increases errors in the estimate. Default 0.
%
% KSUGGESTION is a suggestion for the k:th nearest neighbor 
% that should be used for estimation. The k can't
% be freely set because the estimation algorithm is only defined
% for k > q - 1. Value zero means an accurate (q-dependent) default 
% is used. For accurate results one should choose 
% kSuggestion >= 2 * ceil(q) - 1. Default 0.
%
% THREADS determines the number of threads to use for parallelization.
% To fully take advantage of multiple cores in your machine, set this
% to the number of cores in your machine. Note however that this makes 
% your computer unresponsive to other tasks. When you need responsiveness, 
% spare one core for other work. Default maxNumCompThreads.

% Description: Temporal Renyi entropy estimation
% Detail: Leonenko-Pronzato-Savani nearest neighbor estimator
% Documentation: tim_matlab_matlab.txt

function H = renyi_entropy_kl_t(...
    S, timeWindowRadius, q, epsilon, kSuggestion, threads)

if nargin < 2
    error('Not enough input arguments.');
end

if nargin > 6
    error('Too many input arguments.');
end

if nargout > 1
    error('Too many output arguments.');
end

if nargin < 3
    q = 2;
end

if nargin < 4
    epsilon = 0;
end

if nargin < 5
    kSuggestion = 0;
end

if nargin < 6
    threads = maxNumCompThreads;
end

check_signalset(S);

if size(timeWindowRadius, 1) ~= 1 || ...
   size(timeWindowRadius, 2) ~= 1
    error('TIMEWINDOWRADIUS must be a scalar.');
end

if timeWindowRadius < 0
    error('TIMEWINDOWRADIUS must be non-negative.');
end

if size(q, 1) ~= 1 || ...
   size(q, 2) ~= 1
    error('Q must be a scalar.');
end

if size(epsilon, 1) ~= 1 || ...
   size(epsilon, 2) ~= 1
    error('EPSILON must be a scalar.');
end

if q <= 0
	error('Q must be positive');
end

if epsilon < 0
    error('EPSILON must be non-negative.');
end

if size(kSuggestion, 1) ~= 1 || ...
   size(kSuggestion, 2) ~= 1
    error('KSUGGESTION must be a scalar integer.');
end

if kSuggestion < 0
    error('KSUGGESTION must be non-negative.');
end

if size(threads, 1) ~= 1 || ...
   size(threads, 2) ~= 1
    error('THREADS must be a scalar integer.');
end

if threads < 1
    error('THREADS must be at least 1.');
end

H = tim_matlab('renyi_entropy_lps_t', ...
    S, timeWindowRadius, q, epsilon, kSuggestion, threads);