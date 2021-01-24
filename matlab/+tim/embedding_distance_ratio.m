function distanceRatioSet = embedding_distance_ratio(pointSet, varargin)

import([tim_package, '.*']);

concept_check(nargin, 'inputs', 1);
concept_check(nargout, 'outputs', 0 : 1);

% Optional input arguments.

embeddingFactorSet = 1 : 16;
eval(process_options({'embeddingFactorSet'}, varargin));

d = size(pointSet, 1);
n = size(pointSet, 2);

embeddingFactors = numel(embeddingFactorSet);

distanceRatioSet = zeros(embeddingFactors, n);
for i = 1 : embeddingFactors
    embeddingFactor = embeddingFactorSet(i);
    
    lowerPointSet = tim.delay_embed(pointSet, 1, 1);
    higherPointSet = tim.delay_embed(pointSet, embeddingFactor + 1, 1);

    kdTree = pastelmatlab.PointKdTree(d);
    idSet = kdTree.insert(lowerPointSet);
    kdTree.refine();

    neighborSet = kdTree.search_nearest(idSet, ...
        ones(1, n) * Inf, 1);

    aIdSet = idSet(neighborSet > 0);
    bIdSet = neighborSet(neighborSet > 0);

    lowerDistanceSet = sqrt(sum(...
        (lowerPointSet(:, aIdSet) - lowerPointSet(:, bIdSet)).^2));
    higherDistanceSet = sqrt(sum(...
        (higherPointSet(:, aIdSet) - higherPointSet(:, bIdSet)).^2));
    
    distanceRatioSet(i, :) = higherDistanceSet / lowerDistanceSet;
end
