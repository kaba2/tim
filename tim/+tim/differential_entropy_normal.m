% DIFFERENTIAL_ENTROPY_NORMAL
% Differential entropy of a normal distribution.
%
% H = differential_entropy_normal(d, 'key', value, ...)
%
% where
%
% D is a positive integer which specifies the dimension of
% the distribution.
%
% Optional input arguments in 'key'-value pairs:
%
% DETCOV ('detCov') is a positive real number which specifies 
% the determinant of the covariance matrix of the normal
% distribution.
% Default: 1.

% Description: Differential entropy of a normal distribution
% Documentation: differential_entropy.txt

function H = differential_entropy_normal(d, varargin)

import([tim_package, '.*']);

concept_check(nargin, 'inputs', 1);
concept_check(nargout, 'outputs', 0 : 1);

% Concept checks
pastelsys.concept_check(...
	d, 'integer', ...
	d, 'positive');

% Optional input arguments
detCov = 1;
eval(process_options({'detCov'}, varargin));

C = log(2 * pi) + 1;
H = 0.5 * (log(detCov) + d * C);