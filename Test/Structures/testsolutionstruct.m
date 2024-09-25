classdef testsolutionstruct < matlab.unittest.TestCase
    methods (Test)
        % Test the class constructor and method:
        function constructortest(testCase)
            A = ones(10,3);
            sol = solutionstruct(A);
            testCase.verifyClass(sol,"solutionstruct");
            testCase.verifyEqual(sol.values,A);
        end

        function valuestest(testCase)
            A = ones(10,3);
            sol = solutionstruct(A);
            testCase.verifyEqual(num(sol),A);
        end
    end
end