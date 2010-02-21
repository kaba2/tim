% DIFFERENTIAL_ENTROPY_NK
% A differential entropy estimate from samples
% using Nilsson-Kleijn manifold nearest neighbor
% estimator.
%
% [H, d] = differential_entropy_nk(S, epsilon, threads)
%
% or
%
% H = differential_entropy_nk(S, epsilon, threads)
%
% where
%
% H is the differential entropy estimate of a random distribution
% lieing on a d-dimensional differentiable manifold, where d is
% also estimated from the data and is an integer.
%
% S is an arbitrary-dimensional cell-array whose linearization contains
% q trials of a signal. Each signal is a real (m x n)-matrix that 
% contains n samples of an m-dimensional signal.
%
% EPSILON is the maximum relative error in distance that
% nearest neighbor searching is allowed to result in.
% Higher tolerances result in enhanced performance, but
% increases errors in the estimate. Default 0.
%
% THREADS determines the number of threads to use for parallelization.
% To fully take advantage of multiple cores in your machine, set this
% to the number of cores in your machine. Note however that this makes 
% your computer unresponsive to other tasks. When you need responsiveness, 
% spare one core for other work. Default maxNumCompThreads.

% Description: Differential entropy estimation
% Detail: Nilsson-Kleijn manifold nearest neighbor estimator
% Documentation: tim_matlab_matlab.txt

function [H, d] = differential_entropy_nk(S, epsilon, threads)

if nargin < 1
    error('Not enough input arguments.');
end

if nargin > 3
    error('Too many input arguments.');
end

if nargout > 2
    error('Too many output arguments.');
end

if nargin < 2
    epsilon = 0;
end

if nargin < 3
    threads = maxNumCompThreads;
end

check_signalset(S);

if size(epsilon, 1) ~= 1 || ...
   size(epsilon, 2) ~= 1
    error('EPSILON must be a scalar.');
end

if epsilon < 0
    error('EPSILON must be non-negative.');
end

if size(threads, 1) ~= 1 || ...
   size(threads, 2) ~= 1
    error('THREADS must be a scalar integer.');
end

if threads < 1
    error('THREADS must be at least 1.');
end

[H, dOut] = tim_matlab('differential_entropy_nk', ...
	S, epsilon, threads);

if nargout > 1
    d = dOut;
end
