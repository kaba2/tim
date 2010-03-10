% DIVERGENCE_WKV
% A Kullback-Leibler divergence estimate from samples
% using Wang-Kulkarni-Verdu nearest neighbor estimator.
%
% D = divergence_wkv(X, Y, threads)
%
% where
%
% X and Y are signal sets.
%
% Type 'help tim_matlab' for documentation.

% Description: Kullback-Leibler divergence estimation
% Detail: Wang-Kulkarni-Verdu nearest neighbor estimator
% Documentation: tim_matlab_matlab.txt

function D = divergence_wkv(X, Y, threads)

if nargin < 2
    error('Not enough input arguments.');
end

if nargin > 3
    error('Too many input arguments.');
end

if nargout > 1
    error('Too many output arguments.');
end

if nargin < 3
	threads = maxNumCompThreads;
end

if isnumeric(X)
    D = divergence_wkv({X}, Y, threads);
    return
end

if isnumeric(Y)
    D = divergence_wkv(X, {Y}, threads);
    return
end

check_signalset(X);
check_signalset(Y);

if numel(X) ~= numel(Y)
	error('The number of trials in X and Y do not match.');
end

if size(X{1}, 1) ~= size(Y{1}, 1)
    error('The dimensions of X and Y do not match.');
end

if size(threads, 1) ~= 1 || ...
   size(threads, 2) ~= 1
    error('THREADS must be a scalar integer.');
end

if threads < 1
    error('THREADS must be at least 1.');
end

D = tim_matlab('divergence_wkv', ...
	X, Y, threads);

end
