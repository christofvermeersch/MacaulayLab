classdef testbkk < matlab.unittest.TestCase
    methods (Test)
        % Test the method:
        function thesistest(testCase)
            c = randn(7,1);
            p = [c(1) 3 2; c(2) 1 0; c(3) 0 2; c(4) 0 0];
            q = [c(5) 1 4; c(6) 3 0; c(7) 0 1];
            actual = bkk(systemstruct({p,q}));
            expected = 18;
            testCase.verifyEqual(actual,expected);
        end
    end
end