classdef testshapiro < matlab.unittest.TestCase
    methods (Test)
        % Test the method:
        function oneparametertest(testCase)
            actual = shapiro(randmep(1,1,10,10));
            expected = 10;
            testCase.verifyEqual(actual,expected);
        end

        function twoparametertest(testCase)
            actual = shapiro(randmep(1,2,3,2));
            expected = 3;
            testCase.verifyEqual(actual,expected);
        end

        function threeparametertest(testCase)
            actual = shapiro(randmep(1,3,5,3));
            expected = 10;
            testCase.verifyEqual(actual,expected);
        end

        function thesistest(testCase)
            A00 = [2 6; 4 5; 0 1];
            A10 = [1 0; 0 1; 1 1];
            A01 = [4 2; 0 8; 1 1];
            mep = mepstruct({A00,A01,A10},[0 0; 0 1; 1 0]);
            actual = shapiro(mep);
            expected = 3;
            testCase.verifyEqual(actual,expected);
        end
    end
end
