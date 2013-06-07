% AUTO_MI_PLOT
% Draws a plot of mutual-information with lagged versions of itself.
%
% auto_mi_plot(pointSet)
%
% Optional arguments
% ------------------
%
% MAXLAG ('maxLag') is a non-negative integer which specifies the
% maximum lag for which to plot the auto mutual information.
% Default: 250
%
% LAGS ('lags') is a positive integer which specifies the number
% of lags to evaluate in the [0, maxLag] range. The lags will be
% positioned uniformly on this range. If LAGS is greater than
% MAXLAG + 1, it will be set to MAXLAG + 1.
% Default: 40

function auto_mi_plot(pointSet, varargin)

% Optional input arguments.
maxLag = 250;
lags = 40;
eval(tim.process_options(...
    {'maxLag', 'lags'}, ...
    varargin));

if lags > maxLag + 1
	% Because lags are integers, evaluating
	% more lags than maxLag + 1 is useless;
	% otherwise same lags will be evaluated
	% multiple times.
	lags = maxLag + 1;
end

% Compute the lags at which to evaluate
% auto mutual information.
lagSet = round(linspace(0, maxLag, lags));

% Compute the auto mutual information for the
% given lags.
miSet = tim.mutual_information(...
    pointSet, pointSet, ...
    'yLag', lagSet)';

% Plot the auto mutual information against the lags.
hold on;
plot(lagSet, miSet);
[ignore, peakSet] = findpeaks(-miSet, 'npeaks', 1);
if ~isempty(peakSet)
	% If there are local minima, add a marker for the
	% first local minimum.	
	plot(lagSet(peakSet), miSet(peakSet), 'k^', 'markerfacecolor', 'r');
	title(['Auto mutual information (first local min = ', num2str(lagSet(peakSet(1))), ')']);
else
	title('Auto mutual information');
end
xlabel('Lag');
ylabel('Mutual information');
hold off;
