% SPATIAL_PLOT
% Visualizes a time-series.
%
% spatial_plot(pointSet, 'key', value, ...)
%
% POINTSET is a real (d x n)-array where each column contains
% a d-dimensional point.
%
% If POINTSET is more than three-dimensional, then only
% the first three dimensions are visualized for the point-set.
%
% Optional arguments
% ------------------
%
% SUBSET ('subset'): An array whose linearization contains indices 
% to POINTSET. Only the points in this subset will be visualized.
% Line segments will be drawn only between sub-sequent indices.
% Default: 1 : n
%
% DRAW_POINTS ('draw_points'): A boolean specifying whether the
% the trajectory should be drawn as colored points.
% Default: false
%
% DRAW_LINES ('draw_lines'): A boolean specifying whether the
% the trajectory should be drawn as line segments.
% Default: true

function spatial_plot(pointSet, varargin)

[d, n] = size(pointSet);

% Optional input arguments.
subset = 1 : n;
draw_points = false;
draw_lines = true;
eval(tim.process_options(...
    {'subset', 'draw_points', 'draw_lines'}, ...
    varargin));

pointSize = 5;
colorMap = cool(n);
actualSet = pointSet(:, subset(:));
colorMap = colorMap(subset(:), :);
m = numel(subset(:));

if d == 1
    if draw_points
        scatter(actualSet, zeros(1, numel(subset)), ...
            pointSize, colorMap, 'o','filled');
    end
elseif d == 2
    hold on;
    if draw_lines
        for i = 1 : (m - 1)
            if subset(i + 1) == subset(i) + 1
                plot(...
                    [actualSet(1, i); actualSet(1, i + 1)], ...
                    [actualSet(2, i); actualSet(2, i + 1)], ...
                    'b');
            end
        end
    end
    if draw_points
        scatter(actualSet(1, :), actualSet(2, :), ...
            pointSize, colorMap, 'o', 'filled');
    end
    hold off;
elseif d >= 3
    hold on;
    if draw_lines
        for i = 1 : (m - 1)
            if subset(i + 1) == subset(i) + 1
                plot3(...
                    [actualSet(1, i); actualSet(1, i + 1)], ...
                    [actualSet(2, i); actualSet(2, i + 1)], ...
                    [actualSet(3, i); actualSet(3, i + 1)], 'b');
            end
        end
    end
    if draw_points
        scatter3(actualSet(1, :), actualSet(2, :), actualSet(3, :), ...
            pointSize, colorMap, 'o', 'filled');
    end
    view(3);
    hold off;
end
