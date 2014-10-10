% AUTO_CORRELATION_PLOT
% Draws a plot of auto-correlation with lagged versions of itself.
%
% hPlot = auto_correlation_plot(pointSet)
%
% where
%
% HPLOT is the handle returned by plot().

function auto_correlation_plot(pointSet)

maxLag = 200;

acSet = xcorr(...
    pointSet, pointSet, ...
    maxLag, 'unbiased');
acSet = acSet((numel(acSet) - 1) / 2 + 2 : end);

lagSet = 1 : maxLag;

hold on;
hPlot = plot(lagSet, acSet);
[ignore, peakSet] = findpeaks(-acSet, 'npeaks', 1);
plot(lagSet(peakSet), acSet(peakSet), 'k^', 'markerfacecolor', 'r');
title('Auto correlation');
xlabel('Lag');
ylabel('Correlation');
hold off;
