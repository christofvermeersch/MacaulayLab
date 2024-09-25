classdef testgrnlex < matlab.unittest.TestCase
    methods (Test)
        % Test the method:
        function simpletest(testCase)
            input = [0 0 0; 1 0 0; 0 1 0; 0 0 1; 2 0 0; 1 1 0; 1 0 1; ...
                0 2 0; 0 1 1; 0 0 2];
            actual = grnlex(input);
            expected = [1; 2; 3; 4; 5; 6; 7; 8; 9; 10];
            testCase.verifyEqual(actual,expected);
        end
    end
end