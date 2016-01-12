% DIFFERENTIAL_ENTROPY_UNIFORM
% Differential entropy of a uniform distribution.
%
% H = differential_entropy_uniform(supportVolume, 'key', value, ...)
%
% where
%
% SUPPORTVOLUME ('supportVolume') is a positive real number
% which specifies the volume of the support of the uniform
% distribution.
% Default: 1

% Description: Differential entropy of a uniform distribution
% Documentation: differential_entropy.txt

function H = differential_entropy_uniform(supportVolume, varargin)

import([tim_package, '.*']);

concept_check(nargin, 'inputs', 1);
concept_check(nargout, 'outputs', 0 : 1);

pastelsys.concept_check(...
	supportVolume, 'real', ...
	supportVolume, 'positive');

H = tim_matlab('differential_entropy_uniform', ...
	supportVolume);
