% DENSITY_PLOT
% Draws an empirical probability density plot.
%
% [hPlot, pdfSet] = density_plot(pointSet, edgeSet)
%
% where
%
% POINTSET is a real array whose linearization contains samples
% drawn from a univariate real random variable.
%
% EDGESET is a real array whose linearization specifies the edges
% which specify the bins to use. The k:th bin contains the samples 
% in the range [EDGESET(k), edgeSet(k + 1)).
%
% HPLOT is the handle returned by bar().
%
% PDFSET is the empirical probability density as computed
% by tim.empirical_density(pointSet, edgeSet).

% Description: Draws an empirical probability density plot

function [hPlot, pdfSet] = density_plot(pointSet, edgeSet, varargin)

import([tim_package, '.*']);

concept_check(nargin, 'inputs', 2);
concept_check(nargout, 'outputs', 0 : 2);

% Linearize the edge-set.
edgeSet = edgeSet(:)';

% Compute the empirical probability density.
pdfSet = tim.empirical_density(pointSet, edgeSet);

% Draw the empirical probability density.
hPlot = bar(edgeSet, [pdfSet, 1], 'histc');
shading('flat');
ylabel('Probability density');

% Compute the width of the support.
width = edgeSet(end) - edgeSet(1);

% Compute the height of the window such that it
% provides a view with 1.5 area; this is 1.5 times
% probability 1.
height = (1 / width) * 1.5;

% Set the window. The first and last edges provide
% the support of the probability density, and therefore
% are natural limits for the x-axis. The probability
% density is non-negative, and therefore 0 is a good
% lower limit on the y-axis. The upper y-limit we
% chose above.
axis([edgeSet(1), edgeSet(end), 0, height]);

