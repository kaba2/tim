% DIFFERENTIAL_ENTROPY_NK
% A differential entropy estimate from samples
% using Nilsson-Kleijn manifold nearest neighbor
% estimator.
%
% [H, d] = differential_entropy_nk(S, threads)
%
% or
%
% H = differential_entropy_nk(S, threads)
%
% where
%
% S is a signal set.
%
% H is the differential entropy estimate of a random distribution
% lieing on a d-dimensional differentiable manifold, where d is
% also estimated from the data and is an integer.
%
% Type 'help tim' for more documentation.

% Description: Differential entropy estimation
% Detail: Nilsson-Kleijn manifold nearest neighbor estimator
% Documentation: tim_matlab_matlab.txt

function [H, d] = differential_entropy_nk(S, threads)

if nargin < 1
    error('Not enough input arguments.');
end

if nargin > 2
    error('Too many input arguments.');
end

if nargout > 2
    error('Too many output arguments.');
end

if nargin < 2
    threads = maxNumCompThreads;
end

if isnumeric(S)
    H = differential_entropy_nk({S}, threads);
    return
end


check(S, 'signalSet');
check(threads, 'threads');

[H, dOut] = tim_matlab('differential_entropy_nk', ...
	S, threads);

if nargout > 1
    d = dOut;
end
