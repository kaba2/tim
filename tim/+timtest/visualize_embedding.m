function visualize_embedding(pointSet)

figure;
plot3(pointSet(1, :), pointSet(2, :), pointSet(3, :));
title('Original point-set')

figure;
plot(pointSet(1, :));
title('X-coordinates')

% Find the embedding delay
% ------------------------

% A good embedding lag is given by the first minimum of the
% auto mutual information. If the signal is sampled densely enough, 
% temporally close samples are related simply because of continuity. 
% Increasing the lag decreases the continuity-caused dependency
% between delayed samples. On the other hand, at some point the
% dependency starts to rise again, since a chaotic system usually
% has cyclic-like properties. Overall, the dependency will lower,
% since temporally distant samples have less and less dependency to
% each other, because of the chaotic pseudo-randomness. This is
% why one should choose the first minimum, rather than subsequent
% minima.

maxLag = 200;
miLags = 100;

miLagSet = round(linspace(0, maxLag, miLags));

miSet = tim.mutual_information(...
    pointSet(1, :), pointSet(1, :), ...
    'yLag', miLagSet);

figure;
plot(miLagSet, miSet);
hold on;
[ignore, miPeakSet] = findpeaks(-miSet, 'npeaks', 1);
plot(miLagSet(miPeakSet), miSet(miPeakSet), 'k^', 'markerfacecolor', 'r');
title('Auto mutual information of the X-coordinates');
xlabel('Lag');
ylabel('Mutual information');
hold off;

% Auto correlation is much worse in predicting a good embedding lag.
% The problem is that auto correlation is only sensitive to linear
% dependencies in the data.

acSet = xcorr(...
    pointSet(1, :), pointSet(1, :), ...
    maxLag, 'unbiased');
acSet = acSet((numel(acSet) - 1) / 2 + 1 : end);

acLagSet = 0 : maxLag;

figure;
plot(acLagSet, acSet);
hold on;
[ignore, acPeakSet] = findpeaks(-acSet, 'npeaks', 1);
plot(acLagSet(acPeakSet), acSet(acPeakSet), 'k^', 'markerfacecolor', 'r');
title('Auto correlation of the X-coordinates');
xlabel('Lag');
ylabel('Correlation');
hold off;

% Create the delay-embeddings
% ---------------------------

bestLag = miLagSet(miPeakSet(1));

% 2D delay-embedding.
embeddedSet = tim.delay_embed(pointSet(1, :), 2, bestLag);

lagText = ['(lag ', num2str(bestLag), ')'];

figure;
plot(embeddedSet(1, :), embeddedSet(2, :));
title(['2D delay-embedding of the X-coordinates ', ...
    lagText]);

% 3D delay-embedding.
embeddedSet = tim.delay_embed(pointSet(1, :), 3, bestLag);

figure;
plot3(embeddedSet(1, :), embeddedSet(2, :), embeddedSet(3, :));
title(['3D delay-embedding of the X-coordinates ', ...
    lagText]);


