classdef testshift < matlab.unittest.TestCase
    properties (TestParameter)
        n = {1,2,3};
        d = {1,3,5};
        basis = {@monomial,@chebyshev};
    end
        
    methods(Test)
        function blockedtest(testCase)
            A = randi(5,3,3);
            B = randi(5,3,3);
            actual = shift(A,B);
            expected = monomial(A,B);
            testCase.verifyEqual(actual,expected);
        end
    end

    methods (Test, ParameterCombination = "exhaustive")
        function defaulttest(testCase,n,d)
            A = randi(d,1,n);
            B = randi(d,1,n);
            actual = shift(A,B);
            expected = monomial(A,B);
            testCase.verifyEqual(actual,expected);
        end

        function optiontest(testCase,n,d,basis)
            A = randi(d,1,n);
            B = randi(d,1,n);
            actual = shift(A,B,basis);
            expected = basis(A,B);
            testCase.verifyEqual(actual,expected);
        end
    end
end
