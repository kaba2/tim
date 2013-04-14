function distanceMatrix2 = distance_matrix(pointSet, varargin)
    
import([tim_package, '.*']);

concept_check(nargin, 'inputs', 1);
concept_check(nargout, 'outputs', 0 : 1);

% Optional input arguments.
eval(process_options({}, varargin));

n = size(pointSet, 2);

distanceMatrix2 = zeros(n, n);
for i = 1 : n
    distanceMatrix2(i, (i + 1) : n) = ...
        sum((pointSet(:, i + 1 : n) - pointSet(:, i) * ones(1, n - i)).^2);
    % The distance matrix is symmetric; copy the results.
    distanceMatrix2((i + 1) : n, i) = distanceMatrix2(i, (i + 1) : n);
end


