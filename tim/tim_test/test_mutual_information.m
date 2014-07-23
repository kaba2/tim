classdef test_mutual_information < matlab.unittest.TestCase
    
    properties (TestParameter)
        xd = struct('one', 1, 'two', 2, 'four', 4, 'ten', 10);
        yd = struct('one', 1, 'two', 2, 'four', 4, 'ten', 10);
        k = struct('one', 1, 'two', 2, 'four', 4, 'ten', 10);
    end
    
    methods (Test)       
        function testIndependent(testCase, xd, yd, k)
            n = 10000;
            X = randn(xd, n);
            Y = randn(yd, n);

            mi = tim.mutual_information(...
                X, Y, 'k', k);

            testCase.verifyEqual(mi, 0, 'AbsTol', 0.05);
        end                
    end
    
end
