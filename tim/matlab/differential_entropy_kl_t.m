% DIFFERENTIAL_ENTROPY_KL_T
% A temporal differential entropy estimate from samples
% using Kozachenko-Leonenko nearest neighbor estimator.
%
% H = differential_entropy_kl_t(
%     S, timeWindowRadius, k, filter, threads)
%
% where
%
% S is a signal set.
%
% Type 'help tim' for more documentation.

% Description: Temporal differential entropy estimation
% Detail: Kozachenko-Leonenko nearest neighbor estimator
% Documentation: tim_matlab_matlab.txt

function H = differential_entropy_kl_t(...
    S, timeWindowRadius, k, filter, threads)

if nargin < 2
    error('Not enough input arguments.');
end

if nargin > 5
    error('Too many input arguments.');
end

if nargout > 1
    error('Too many output arguments.');
end

if nargin < 3
    k = 1;
end

if nargin < 4
    filter = [1];
end

if nargin < 5
    threads = maxNumCompThreads;
end

if isnumeric(S)
    H = differential_entropy_kl_t(...
        {S}, timeWindowRadius, k, filter, threads);
    return
end

check(S, 'signalSet');
check(timeWindowRadius, 'timeWindowRadius');
check(k, 'k');
check(filter, 'filter');
check(threads, 'threads');

H = tim_matlab('differential_entropy_kl_t', ...
    S, timeWindowRadius, k, filter, threads);
