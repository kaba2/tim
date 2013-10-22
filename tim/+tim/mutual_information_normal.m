% MUTUAL_INFORMATION_NORMAL
% Mutual information between normal distributions.
%
% I = mutual_information_normal(xDetCov, yDetCov, xyDetCov, ...
%                               'key', value, ...)
%
% XDETCOV is a positive real number which specifies the determinant
% of the covariance matrix of the random variable X. 
%
% YDETCOV is a positive real number which specifies the determinant
% of the covariance matrix of the random variable Y.
%
% XYDETCOV is a positive real number which specifies the determinant
% of the covariance matrix of the random variable (X, Y).

% Description: Mutual information between normal distributions.
% Documentation: mutual_information.txt

function I = mutual_information_normal(...
    xDetCov, yDetCov, xyDetCov, varargin)

import([tim_package, '.*']);

pastelsys.concept_check(...
	xDetCov, 'real', ...
	xDetCov, 'positive', ...
	yDetCov, 'real', ...
	yDetCov, 'positive', ...
	xyDetCov, 'real', ...
	xyDetCov, 'positive')

% Optional input arguments
eval(process_options({}, varargin));

% This is the mutual information evaluated by differential
% entropies:
%
%    I(X, Y) = H(X) + H(Y) - H(X, Y)
%
% See differential_entropy_normal.m for the analytical formula
% for the differential entropy of the normal distribution.
% As a check, if the covariance matrix of (X, Y) is block-diagonal,
% with cov((X, Y)) = diag(cov(X), cov(Y)), that is, X and Y are
% independent, then det(cov((X, Y))) = det(cov(X)) det(cov(Y)), and
% the mutual information is zero.

I = 0.5 * (log(xDetCov) + log(yDetCov) - log(xyDetCov));
