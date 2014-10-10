% AUTO_MI_PLOT
% Draws a plot of mutual-information with lagged versions of itself.
%
% hPlot = auto_mi_plot(pointSet)
% hPlot = auto_mi_plot(pointSet, 'key', value, ...)
%
% where
%
% POINTSET is a real (d x n)-matrix whose columns contain
% d-dimensional points.
%
% HPLOT is the handle returned by plot().
%
% Optional arguments
% ------------------
%
% LAGSET ('lagSet') is an integer array whose linearization 
% contains p lags. Default: [0 : 63]

function hPlot = auto_mi_plot(pointSet, varargin)

import([tim_package, '.*']);

concept_check(nargin, 'inputs', 1);
concept_check(nargout, 'outputs', 0 : 1);

% Optional input arguments.
lagSet = 0 : 63;
eval(tim.process_options(...
    {'lagSet'}, ...
    varargin));

% Compute the auto mutual information.
miSet = tim.auto_mi(pointSet, 'lagSet', lagSet);

% Plot the auto mutual information against the lags.
hold on;
hPlot = plot(lagSet, miSet);
[ignore, peakSet] = findpeaks(-miSet, 'npeaks', 1);
if ~isempty(peakSet)
	% If there are local minima, add a marker for the
	% first local minimum.	
	plot(lagSet(peakSet), miSet(peakSet), 'k^', 'markerfacecolor', 'r');
	title(['Auto mutual information (first local min = ', ...
        num2str(lagSet(peakSet(1))), ')']);
else
	title('Auto mutual information');
end
xlabel('Lag');
ylabel('Mutual information');
axis([lagSet(1), lagSet(end), 0, ceil(max(miSet))]);
hold off;
