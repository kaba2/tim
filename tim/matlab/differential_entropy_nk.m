% DIFFERENTIAL_ENTROPY_NK
% A differential entropy estimate from samples
% using Nilsson-Kleijn manifold nearest neighbor
% estimator.
%
% [H, d] = differential_entropy_nk(S)
%
% or
%
% H = differential_entropy_nk(S)
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

function [H, d] = differential_entropy_nk(S)

check(nargin, 'inputs', 1);
check(nargout, 'outputs', 0 : 2);

if isnumeric(S)
    S = {S};
end


check(S, 'signalSet');

[H, dOut] = tim_matlab(...
    'differential_entropy_nk', ...
	S);

if nargout > 1
    d = dOut;
end
