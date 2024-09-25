classdef testgrvlex < matlab.unittest.TestCase
    methods (Test)
        % Test the method:
        function simpletest(testCase)
            input = [0 0 0; 1 0 0; 0 1 0; 0 0 1; 2 0 0; 1 1 0; 1 0 1; ...
                0 2 0; 0 1 1; 0 0 2];
            actual = grvlex(input);
            expected = [1; 2; 3; 4; 5; 6; 8; 7; 9; 10];
            testCase.verifyEqual(actual,expected);
        end
    end
end