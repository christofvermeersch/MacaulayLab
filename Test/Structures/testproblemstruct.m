classdef testproblemstruct < matlab.unittest.TestCase
    methods (Test)
        % Test the class constructor:
        function constructortest(testCase)
            coef = {randn(4,3,2);randn(4,3,2)};
            supp = {[0 0; 1 0; 0 1; 2 0];[0 0; 2 2; 0 1; 0 5]};
            problem = problemstruct(coef,supp);
            testCase.verifyClass(problem,"problemstruct");
            testCase.verifyEqual(problem.coef,coef);
            testCase.verifyEqual(problem.supp,supp);
            testCase.verifyEqual(problem.n,2);
            testCase.verifyEqual(problem.s,2);
            testCase.verifyEqual(problem.k,3);
            testCase.verifyEqual(problem.l,2);
            testCase.verifyEqual(problem.di,{2;5});
            testCase.verifyEqual(problem.dmax,5);
            testCase.verifyEqual(problem.nnze,{24;24});
            testCase.verifyEqual(problem.basis,"unknown");
        end

        % Test the class functions:
        function propertiestest(testCase)
            coef = {[1; 2];[3; 4]};
            supp = {[5 0; 2 2];[4 2; 1 1]};
            problem = problemstruct(coef,supp);
            [s,dmax,n,di,nnze,k,l] = properties(problem);
            expected = {2,6,2,{5;6},{2;2},1,1};
            testCase.verifyEqual({s,dmax,n,di,nnze,k,l},expected);
        end
    end
end