% VISUALIZE_EMBEDDING
% Visualizes a delay-embedding reconstruction.
%
% visualize_embedding(pointSet, 'key', value, ...)
%
% POINTSET is a real (d x n)-array where each column contains
% a d-dimensional point.
%
% Optional arguments
% ------------------
%
% TDELTA ('tDelta'): An integer specifying the embedding delay.
% Default: First minimum of auto mutual information.
% 
% AXIS ('axis'): An integer specifying the coordinate axis to
% use to reconstruct the signal by delay-embedding.
% Default: 0
%
% If POINTSET is more than three-dimensional, then only
% the first three dimensions are visualized for the point-set.
% Delay-embeddings are visualized for 2D and 3D.

function visualize_embedding(pointSet, varargin)

% Optional input arguments.
tDelta = {};
axis = 1;
eval(tim.process_options(...
    {'tDelta', 'axis'}, ...
    varargin));

n = size(pointSet, 2);
d = size(pointSet, 1);
pointSize = 5;
colorMap = cool(n);

figure;
if d == 2
    scatter(pointSet(1, :), pointSet(2, :), ...
        pointSize, colorMap, 'o', 'filled');
elseif d >= 3
    scatter3(pointSet(1, :), pointSet(2, :), pointSet(3, :), ...
        pointSize, colorMap, 'o', 'filled');
end
title('Original point-set')

figure;
plot(pointSet(axis, :));
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

miLagSet = round(linspace(1, maxLag, miLags));

miSet = tim.mutual_information(...
    pointSet(axis, :), pointSet(axis, :), ...
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

if iscell(tDelta)
    tDelta = miLagSet(miPeakSet(1));
end

% Auto correlation is much worse in predicting a good embedding lag.
% The problem is that auto correlation is only sensitive to linear
% dependencies in the data.

acSet = xcorr(...
    pointSet(axis, :), pointSet(axis, :), ...
    maxLag, 'unbiased');
acSet = acSet((numel(acSet) - 1) / 2 + 2 : end);

acLagSet = 1 : maxLag;

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

% 2D delay-embedding.
embeddedSet = tim.delay_embed(pointSet(axis, :), 2, tDelta);

lagText = ['(lag ', num2str(tDelta), ')'];

figure;
scatter(embeddedSet(1, :), embeddedSet(2, :), ...
    pointSize, colorMap, 'o', 'filled');
title(['2D delay-embedding of the X-coordinates ', ...
    lagText]);

% 3D delay-embedding.
embeddedSet = tim.delay_embed(pointSet(axis, :), 3, tDelta);

figure;
scatter3(embeddedSet(1, :), embeddedSet(2, :), embeddedSet(3, :), ...
    pointSize, colorMap, 'o', 'filled');
title(['3D delay-embedding of the X-coordinates ', ...
    lagText]);


