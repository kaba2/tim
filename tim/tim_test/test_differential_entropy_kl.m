classdef test_differential_entropy_kl < matlab.unittest.TestCase
    
    properties (TestParameter)
        d = struct('one', 1, 'two', 2, 'four', 4, 'ten', 10);
        k = struct('one', 1, 'two', 2, 'four', 4, 'ten', 10);
    end
    
    methods (Test)       
        function testAgainstNormal(testCase, d, k)
            n = 10000;
            pointSet = randn(d, n);

            de = tim.differential_entropy_kl(...
                pointSet, 'k', k);

            % The estimator is very accurate for the standard
            % normal distribution.
            correct = tim.differential_entropy_normal(d);
            testCase.verifyEqual(de, correct, 'RelTol', 0.03);
        end                

        function testAgainstUniform(testCase, d, k)
            n = 10000;
            pointSet = rand(d, n);

            de = tim.differential_entropy_kl(...
                pointSet, 'k', k);

            % The estimator is quite inaccurate for the
            % uniform distribution. This is probably because
            % the distribution is discontinuous, and thus
            % the assumption of a local uniform distribution
            % does not hold at the boundaries.
            correct = tim.differential_entropy_uniform(d);
            testCase.verifyEqual(de, correct, 'AbsTol', 1.5);
        end                
    end
    
end
