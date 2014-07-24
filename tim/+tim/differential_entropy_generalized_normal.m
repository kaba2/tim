% DIFFERENTIAL_ENTROPY_GENERALIZED_NORMAL
% Differential entropy of a generalized normal distribution.
%
% H = differential_entropy_generalized_normal(dimension, 'key', value, ...)
%
% where
%
% DIMENSION is a positive integer which specifies the 
% dimension of the distribution.
%
% Optional input arguments in 'key'-value pairs:
%
% SHAPE ('shape') is a positive real number which specifies the 
% shape-parameter of the distribution.
% Default: sqrt(2)
%
% SCALE ('scale') is a positive real number which specifies the 
% scale-parameter of the distribution.
% Default: 2

% Description: Differential entropy of a generalized normal distribution
% Documentation: differential_entropy.txt

function H = differential_entropy_generalized_normal(dimension, varargin)

import([tim_package, '.*']);

concept_check(nargin, 'inputs', 1);
concept_check(nargout, 'outputs', 0 : 1);

% Concept checks
pastelsys.concept_check(...
	dimension, 'integer', ...
	dimension, 'positive');

% Optional input arguments
shape = sqrt(2);
scale = 2;
eval(process_options({'shape', 'scale'}, varargin));

pastelsys.concept_check(...
	shape, 'real', ...
	shape, 'positive', ...
	scale, 'real', ...
	scale, 'positive');

H = tim_matlab('differential_entropy_generalized_normal', ...
	dimension, shape, scale);
