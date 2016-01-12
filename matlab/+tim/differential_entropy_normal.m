% DIFFERENTIAL_ENTROPY_NORMAL
% Differential entropy of a normal distribution.
%
% H = differential_entropy_normal(dimension, 'key', value, ...)
%
% where
%
% DIMENSION is a positive integer which specifies the 
% dimension of the distribution.
%
% Optional input arguments in 'key'-value pairs:
%
% DETCOV ('detCov') is a positive real number which specifies 
% the determinant of the covariance matrix of the normal
% distribution.
% Default: 1.

% Description: Differential entropy of a normal distribution
% Documentation: differential_entropy.txt

function H = differential_entropy_normal(dimension, varargin)

import([tim_package, '.*']);

concept_check(nargin, 'inputs', 1);
concept_check(nargout, 'outputs', 0 : 1);

% Concept checks
pastelsys.concept_check(...
	dimension, 'integer', ...
	dimension, 'positive');

% Optional input arguments
detCov = 1;
eval(process_options({'detCov'}, varargin));

H = tim_matlab('differential_entropy_normal', ...
	dimension, detCov);
