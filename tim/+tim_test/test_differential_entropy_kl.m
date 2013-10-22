classdef test_differential_entropy_kl < handle
    methods (Static)       
        function input = generate()
            pointSet = randn(10, 10000);
            kSet = 1 : 16;
            input = struct(...
                'pointSet', pointSet, ...
                'kSet', kSet);
        end
        
        function output = test_1_0_0(input)
            hSet = zeros(1, numel(input.kSet));
            for k = input.kSet
                hSet(k) = differential_entropy_kl(...
                    input.pointSet, k);
            end
            output = struct('hSet', hSet);
        end
        
        function output = test_1_2_0(input)
            hSet = zeros(1, numel(input.kSet));
            for k = input.kSet
                hSet(k) = tim.differential_entropy_kl(...
                    input.pointSet, 'k', k);
            end
            output = struct('hSet', hSet);
        end
        
        function valid = validate(input, output)
            d = size(input.pointSet, 1);
            correct = tim.differential_entropy_normal(d);
            relativeSet = abs(output.hSet - correct) / correct;
            valid = all(relativeSet < 0.01);
        end
        
        function equal = crossValidate(input, aOutput, bOutput)
            equal = all(aOutput.hSet == bOutput.hSet);
        end
                
    end
end