classdef testmonomial < matlab.unittest.TestCase
    properties (TestParameter)
        n = {1,2,3};
        nMonomials = {1,2,3};
        nShifts = {1,2,3};
    end
        
    methods (Test)
        % Test the method:
        function shifttest(testCase)
            A = [1 2 3; 4 0 0];
            B = [1 1 1; 2 0 2];
            actual = monomial(A,B);
            expected = [1 2 3 4; 1 3 2 5; 1 5 1 1; 1 6 0 2];
            testCase.verifyEqual(actual,expected);
        end

        function recursiontest(testCase)
            [~,alpha,beta,phi0,phi1] = monomial(ones(1,3),ones(1,3));
            expectedAlpha = [1 1 1 1];
            expectedBeta = [0 0 0 0];
            expectedPhi0 = [1 0 0 0];
            expectedPhi1 = [1 1 1 1];
            testCase.verifyEqual(alpha,expectedAlpha);
            testCase.verifyEqual(beta,expectedBeta);
            testCase.verifyEqual(phi0,expectedPhi0);
            testCase.verifyEqual(phi1,expectedPhi1);
        end
    end

    methods (Test, ParameterCombination = "exhaustive")
        function sizetest(testCase,n,nMonomials,nShifts)
            A = randi(5,nShifts,n);
            B = randi(5,nMonomials,n);
            actual = monomial(A,B);
            expectedSize = [nShifts*nMonomials, n+1];
            testCase.verifySize(actual,expectedSize);
        end
    end
end
