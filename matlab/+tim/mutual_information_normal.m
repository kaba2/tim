% MUTUAL_INFORMATION_NORMAL
% Mutual information between marginals of a normal distribution.
%
% I = mutual_information_normal(productDetCov, jointDetCov, ...
%                               'key', value, ...)
%
% PRODUCTDETCOV is a positive real number which specifies the 
% product of the determinants of the covariance matrices of the 
% normally distributed random variables X_i, where i in [1, p]. 
%
% JOINTDETCOV is a positive real number which specifies the determinant
% of the covariance matrix of the normally distributed random variable 
% X = (X_1, ..., X_p).
%
% When p > 2, the mutual information computed here is called 
% total correlation. The most common case is to have p = 2.

% Description: Mutual information between marginals of a normal distribution.
% Documentation: mutual_information.txt

function I = mutual_information_normal(...
    productDetCov, jointDetCov, varargin)

import([tim_package, '.*']);

concept_check(nargin, 'inputs', 2);
concept_check(nargout, 'outputs', 0 : 1);

pastelmatlab.concept_check(...
	productDetCov, 'real', ...
	productDetCov, 'positive', ...
	jointDetCov, 'real', ...
	jointDetCov, 'positive');

I = tim_matlab('mutual_information_normal', ...
	productDetCov, jointDetCov);
