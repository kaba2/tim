% RENYI_ENTROPY_LPS_T
% A temporal Renyi entropy estimate from samples
% using Leonenko-Pronzato-Savani nearest neighbor estimator.
%
% H = renyi_entropy_lps_t(
%     S, timeWindowRadius, q, kSuggestion, filter, threads)
%
% where
%
% S is a signal set.
%
% Q is the power in the definition Renyi entropy.
% If Q = 1, differential_entropy_kl_t() is used to
% compute the result instead. 
% If Q < 1, there are huge errors in the estimation.
% Default 2.
%
% KSUGGESTION is a suggestion for the k:th nearest neighbor 
% that should be used for estimation. The k can't
% be freely set because the estimation algorithm is only defined
% for k > q - 1. Value zero means an accurate (q-dependent) default 
% is used. For accurate results one should choose 
% kSuggestion >= 2 * ceil(q) - 1. Default 0.
%
% Type 'help tim' for more documentation.

% Description: Temporal Renyi entropy estimation
% Detail: Leonenko-Pronzato-Savani nearest neighbor estimator
% Documentation: tim_matlab_matlab.txt

function H = renyi_entropy_lps_t(...
    S, timeWindowRadius, q, kSuggestion, filter, threads)

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
    kSuggestion = 0;
end

if nargin < 5
    filter = [1];
end

if nargin < 6
    threads = maxNumCompThreads;
end

if isnumeric(S)
    H = renyi_entropy_lps_t({S}, timeWindowRadius, ...
        q, kSuggestion, filter, threads);
    return
end

check(S, 'signalSet');
check(timeWindowRadius, 'timeWindowRadius');

if size(q, 1) ~= 1 || ...
   size(q, 2) ~= 1
    error('Q must be a scalar.');
end

if q <= 0
	error('Q must be positive');
end

if size(kSuggestion, 1) ~= 1 || ...
   size(kSuggestion, 2) ~= 1
    error('KSUGGESTION must be a scalar integer.');
end

if kSuggestion < 0
    error('KSUGGESTION must be non-negative.');
end

check(filter, 'filter');
check(threads, 'threads');

H = tim_matlab('renyi_entropy_lps_t', ...
    S, timeWindowRadius, q, kSuggestion, filter, threads);
