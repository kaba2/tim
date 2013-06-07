% AUTO_CORRELATION_PLOT
% Draws a plot of auto-correlation with lagged versions of itself.
%
% auto_correlation_plot(pointSet)

function auto_correlation_plot(pointSet)

maxLag = 200;

acSet = xcorr(...
    pointSet, pointSet, ...
    maxLag, 'unbiased');
acSet = acSet((numel(acSet) - 1) / 2 + 2 : end);

lagSet = 1 : maxLag;

hold on;
plot(lagSet, acSet);
[ignore, peakSet] = findpeaks(-acSet, 'npeaks', 1);
plot(lagSet(peakSet), acSet(peakSet), 'k^', 'markerfacecolor', 'r');
title('Auto correlation');
xlabel('Lag');
ylabel('Correlation');
hold off;
