classdef test_differential_entropy_kl < matlab.unittest.TestCase
    
    properties (TestParameter)
        d = struct('one', 1, 'two', 2, 'four', 4, 'ten', 10);
        k = struct('one', 1, 'two', 2, 'four', 4, 'ten', 10);
    end
    
    methods (Test)       
        function testAgainstAnalyticNormal(testCase, d, k)
            n = 10000;
            pointSet = randn(d, n);

            h = tim.differential_entropy_kl(...
                pointSet, 'k', k);

            correct = tim.differential_entropy_normal(d);
            relativeError = abs(h - correct) / correct;
            testCase.verifyTrue(relativeError < 0.02);
        end
                
    end
end