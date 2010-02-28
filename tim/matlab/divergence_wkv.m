% DIVERGENCE_WKV
% A Kullback-Leibler divergence estimate from samples
% using Wang-Kulkarni-Verdu nearest neighbor estimator.
%
% D = divergence_wkv(X, Y, threads)
%
% where
%
% X and Y are arbitrary-dimensional cell-arrays whose linearizations 
% contain q trials of signals X and Y, respectively.
%
% THREADS determines the number of threads to use for parallelization.
% To fully take advantage of multiple cores in your machine, set this
% to the number of cores in your machine. Note however that this makes 
% your computer unresponsive to other tasks. When you need responsiveness, 
% spare one core for other work. Default maxNumCompThreads.
%
% Each signal is a real (m x n)-matrix that contains n samples of an
% m-dimensional signal. The signals contained in X (Y) must all have equal
% dimensionality, but their number of samples may vary. 
% The number of trials in X and Y must be equal.

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
