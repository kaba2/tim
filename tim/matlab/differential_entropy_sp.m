% DIFFERENTIAL_ENTROPY_SP
% A differential entropy estimate from samples
% using Stowell-Plumbley recursive partition estimator.
%
% H = differential_entropy_sp(S)
%
% where
%
% S is a signal set.
%
% Type 'help tim_matlab' for documentation.

% Description: Differential entropy estimation
% Detail: Stowell-Plumbley recursive partition estimator
% Documentation: tim_matlab_matlab.txt

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

if isnumeric(S)
    H = differential_entropy_sp({S});
    return
end

check_signalset(S);

if (size(S{1}, 1) > 3)
	warning('This estimator has bad accuracy for dimensions > 3!');
end

H = tim_matlab('differential_entropy_sp', S);
