% EMPIRICAL_DENSITY
% Computes an empirical probability density by binning.
%
% hPlot = empirical_density(pointSet, edgeSet)
%
% where
%
% POINTSET is a real array whose linearization contains samples
% drawn from a univariate real random variable.
%
% EDGESET is a real array whose linearization specifies the edges
% which specify the bins to use. The k:th bin contains the samples 
% in the range [EDGESET(k), edgeSet(k + 1)).

% Description: Computes an empirical probability density by binning

function pdfSet = empirical_density(pointSet, edgeSet, varargin)

import([tim_package, '.*']);

concept_check(nargin, 'inputs', 2);
concept_check(nargout, 'outputs', 0 : 1);

% Linearize the point-set.
pointSet = pointSet(:)';
% Linearize the edge-set.
edgeSet = edgeSet(:)';

% Compute the widths of the bins.
widthSet = diff(edgeSet);

% Compute the number of points inside each bin.
histogram = histc(pointSet, edgeSet);
% The last bin contains the samples which are equal to edgeSet(end).
% We discard these samples.
histogram(end) = [];

% An empirical probability density is obtained from the histogram by
% * normalizing with the total number of points, which gives a 
% probability mass function, and
% * normalizing with the width of the bins, so that the area
% of the bins sum to 1.
pdfSet = histogram ./ (widthSet * sum(histogram));

