% Predictor

% Description: Non-linear predictor (Matlab)

classdef Predictor < handle
    properties (Access = 'private')
        kdTree_
        correlationTime_
        kNearest_
    end
    methods
        function self = Predictor(dimension)
            if dimension <= 0
                error('dimension must be a positive integer');
            end
            
            import pastelgeometry.PointKdTree;

            self.kdTree_ = PointKdTree(dimension);
            self.correlationTime_ = 1;
            self.kNearest_ = 4;
        end

        function setKNearest(self, kNearest)
            if kNearest <= 0
                error('kNearest must positive.')
            end
            self.kNearest_ = kNearest;
        end

        function kNearest_ = kNearest(self)
            kNearest_ = self.kNearest_;
        end

        function setCorrelationTime(self, correlationTime)
            if correlationTime <= 0
                error('correlationTime must positive.')
            end
            self.correlationTime_ = correlationTime;
        end

        function correlationTime_ = correlationTime(self)
            correlationTime_ = self.correlationTime_;
        end

        function insert(self, point)
            id = self.kdTree_.insert(point);
            self.kdTree_.hide(id);
            if id > self.correlationTime_
                self.kdTree_.show(id - self.correlationTime_)
            end

            if mod(self.kdTree_.points(), 100) == 99
                % Recompute the subdivision for better
                % performance.
                self.kdTree_.merge();
                self.kdTree_.refine();
            end
        end

        function predictedSet = predict(self, point)
            currentSet = self.kdTree_.search_nearest(...
                point, inf, self.kNearest);
            futureSet = currentSet + 1;
            predictedSet = self.kdTree_.as_points(futureSet);
        end
    end
end

