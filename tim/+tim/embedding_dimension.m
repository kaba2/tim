function ed = embedding_dimension(pointSet)

import([tim_package, '.*']);

d = size(pointSet, 1);
n = size(pointSet, 2);

kdTree = pastelgeometry.PointKdTree(d);
kdTree.insert(pointSet);
kdTree.refine();
clear kdTree;
