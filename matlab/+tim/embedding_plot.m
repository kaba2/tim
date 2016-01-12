% EMBEDDING_PLOT
% Visualizes a delay-embedding reconstruction.
%
% embedding_plot(pointSet, 'key', value, ...)
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

% Description: Visualizes a delay-embedding reconstruction

function embedding_plot(pointSet, varargin)

% Optional input arguments.
tDelta = {};
axis = 1;
eval(tim.process_options(...
    {'tDelta', 'axis'}, ...
    varargin));

[d, n] = size(pointSet);
colorMap = cool(n);
pointSize = 5;

figure;
tim.random_projection_plot(pointSet);
title('Original point-set')

figure;
plot(pointSet(axis, :));
title('X-coordinates')

% Find the embedding delay
% ------------------------

figure;
tim.auto_mi_plot(pointSet(axis, :));

figure
tim.auto_correlation_plot(pointSet(axis, :));

% Create the delay-embeddings
% ---------------------------

tDelta = tim.embedding_delay(pointSet(axis, :));

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


