classdef testposition < matlab.unittest.TestCase
    properties (TestParameter)
        n = {1,2,3};
        d = {1,3,5};
        order = {@grevlex,@grlex,@grnlex,@grvlex};
    end
        
    methods(Test)
        function blockedtest(testCase)
            monomial = randi(5,3,3);
            actual = position(monomial);
            expected = grevlex(monomial);
            testCase.verifyEqual(actual,expected);
        end
    end

    methods (Test, ParameterCombination = "exhaustive")
        function defaulttest(testCase,n,d)
            monomial = randi(d,1,n);
            actual = position(monomial);
            expected = grevlex(monomial);
            testCase.verifyEqual(actual,expected);
        end

        function optiontest(testCase,n,d,order)
            monomial = randi(d,1,n);
            actual = position(monomial,order);
            expected = order(monomial);
            testCase.verifyEqual(actual,expected);
        end
    end
end
