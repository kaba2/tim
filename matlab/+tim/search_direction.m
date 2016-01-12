function indexSet = search_direction(pointSet, direction, ...
    minDot, maxDot)

n = size(pointSet, 2);

direction = direction(:) / norm(direction);

size(direction)
size(pointSet, 1)

indexSet = [];
for i = 1 : (n - 1)
    iDirection = pointSet(:, i + 1) - pointSet(:, i);
    metric = dot(direction, iDirection / norm(iDirection));
    if metric >= minDot && metric <= maxDot
        indexSet(end + 1) = i;
    end
end
