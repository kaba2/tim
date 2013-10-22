% RANDOM_NORMAL
% Generates samples from a normal distribution.
%
% X = random_normal(d, n, 'key', value, ...)
%
% where
%
% D is a positive integer which specifies the dimension of
% the distribution.
%
% N is a non-negative integer which specifies the number
% of samples to generate.
%
% Optional input arguments in 'key'-value pairs:
%
% COV ('cov') is a positive (D x D) real matrix which specifies 
% the covariance matrix of the normal distribution.
% Default: eye(D, D).

% Description: Generates samples from a normal distribution.

function X = random_normal(d, n, varargin)

import([tim_package, '.*']);

concept_check(nargin, 'inputs', 2);
concept_check(nargout, 'outputs', 0 : 1);

% Concept checks
pastelsys.concept_check(...
	d, 'integer', ...
	d, 'positive', ...
	n, 'integer', ...
	n, 'non_negative');

% Optional input arguments
cov = eye(d, d);
eval(process_options({'cov'}, varargin));

pastelsys.concept_check(...
	cov, 'real_matrix', ...
	cov, 'square_matrix');

L = chol(cov, 'lower');
X = L * randn(d, n);
