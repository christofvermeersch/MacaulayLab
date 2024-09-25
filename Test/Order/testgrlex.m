classdef testgrlex < matlab.unittest.TestCase
    methods (Test)
        % Test the method:
        function simpletest(testCase)
            input = [0 0 0; 1 0 0; 0 1 0; 0 0 1; 2 0 0; 1 1 0; 1 0 1; ...
                0 2 0; 0 1 1; 0 0 2];
            actual = grlex(input);
            expected = [1; 4; 3; 2; 10; 9; 8; 7; 6; 5];
            testCase.verifyEqual(actual,expected);
        end
    end
end