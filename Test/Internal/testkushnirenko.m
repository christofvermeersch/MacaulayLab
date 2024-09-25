classdef testkushnirenko < matlab.unittest.TestCase
    methods (Test)
        % Test the method:
        function thesistest(testCase)
            p = [23 16 14; -1 22 18; 27 0 0];
            q = [35 16 14; -1 22 18; 9 0 0];
            actual = kushnirenko(systemstruct({p,q}));
            expected = 20;
            testCase.verifyEqual(actual,expected);
        end
    end
end