classdef testcompanion < matlab.unittest.TestCase
    properties (Constant)
        tolerance = 1e-10;
    end

    methods (Test)
        % Test the method:
        function matrixtest(testCase)
            actual = companion({[3 1; 2 2; 1 3]});
            expected = [0 0 0; 1 0 -3; 0 1 -2];
            testCase.verifyEqual(actual,expected);
        end

        function eigtest(testCase)
            P = [-6 0; 1 4; -5 3; 5 2; 5 1];
            difference = sum(sort(eig(companion({P}))) - [-1; 1; 2; 3]);
            testCase.verifyLessThanOrEqual(difference,testCase.tolerance);
        end
        
        function complexmatrixtest(testCase)
            actual = companion({[i 0; 1+i 1; 2 2; 1-i 3; 1 4]});
            expected = [0 0 0 -i; 1 0 0 -1-i; 0 1 0 -2; 0 0 1 -1+i];
            testCase.verifyEqual(actual,expected);
        end

        function complexeigtest(testCase)
            P = [42+42i 2; 7+i 1];
            difference = sum(sort(eig(companion({P}))) ...
                - [0; -2/21 + 1i/14]);
            testCase.verifyLessThanOrEqual(difference,testCase.tolerance);
        end
    end
end
