% DIFFERENTIAL_ENTROPY_SP
% A differential entropy estimate from samples
% using Stowell-Plumbley recursive partition estimator.
%
% H = differential_entropy_sp(S)
%
% where
%
% S is an arbitrary dimensional cell array whose linearization contains
% q trials of a signal. Each signal is a real (m x n)-matrix that 
% contains n samples of an m-dimensional signal.

% Description: Differential entropy estimation
% Detail: Stowell-Plumbley recursive partition estimator
% Documentation: tim_matlab.txt

function H = differential_entropy_sp(S)

% The limit for the dimension is arbitrary, but
% protects for the case when the user accidentally
% passes the transpose of the intended data.
maxDimension = 32;

if nargin < 1
    error('Not enough input arguments.');
end

if nargin > 1
    error('Too many input arguments.');
end

if nargout > 1
    error('Too many output arguments.');
end

if ~iscell(S)
    error('S must be a cell array.');
end

signals = prod(size(S));

if signals < 1
    error('S must contain at least 1 signal.');
end

dimension = size(S{1}, 1);

for i = 1 : signals
    if ~isa(S{i}, 'double')
        error('Some signal of S is not of type double.');
    end

    if size(S{i}, 1) ~= dimension
        error(['The dimensions of the trials do not match.']);
    end
    
    if size(S{i}, 1) > maxDimension
        error(['Some signal of S has dimension greater than ', ...
            int2str(maxDimension), '.']);
    end
end

H = timDifferentialEntropySp(S);
