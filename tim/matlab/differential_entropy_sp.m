% DIFFERENTIAL_ENTROPY_SP
% A differential entropy estimate from samples
% using Stowell-Plumbley recursive partition estimator.
%
% H = differential_entropy_sp(S)
%
% where
%
% S is an arbitrary-dimensional cell-array whose linearization contains
% q trials of a signal. Each signal is a real (m x n)-matrix that 
% contains n samples of an m-dimensional signal.

% Description: Differential entropy estimation
% Detail: Stowell-Plumbley recursive partition estimator
% Documentation: tim_matlab.txt

function H = differential_entropy_sp(S)

if nargin < 1
    error('Not enough input arguments.');
end

if nargin > 1
    error('Too many input arguments.');
end

if nargout > 1
    error('Too many output arguments.');
end

check_signalset(S);

if (size(S{1}, 1) > 3)
	warning('This estimator has bad accuracy for dimensions > 3!');
end

H = tim_matlab('differential_entropy_sp', S);
