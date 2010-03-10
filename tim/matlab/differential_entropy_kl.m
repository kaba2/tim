% DIFFERENTIAL_ENTROPY_KL
% A differential entropy estimate from samples
% using Kozachenko-Leonenko nearest neighbor estimator.
%
% H = differential_entropy_kl(S, k, threads)
%
% where
%
% S is a signal set.
%
% Type 'help tim' for more documentation.

% Description: Differential entropy estimation
% Detail: Kozachenko-Leonenko nearest neighbor estimator
% Documentation: tim_matlab_matlab.txt

function H = differential_entropy_kl(S, k, threads)

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
    threads = maxNumCompThreads;
end

if isnumeric(S)
    H = differential_entropy_kl({S}, k, threads);
    return
end

check_signalset(S);

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

H = tim_matlab('differential_entropy_kl', ...
	S, k, threads);
