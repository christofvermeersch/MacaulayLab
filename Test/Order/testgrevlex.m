classdef testgrevlex < matlab.unittest.TestCase
    methods (Test)
        % Test the method:
        function simpletest(testCase)
            input = [0 0 0; 1 0 0; 0 1 0; 0 0 1; 2 0 0; 1 1 0; 1 0 1; ...
                0 2 0; 0 1 1; 0 0 2];
            actual = grevlex(input);
            expected = [1; 4; 3; 2; 10; 9; 7; 8; 6; 5];
            testCase.verifyEqual(actual,expected);
        end
    end
end