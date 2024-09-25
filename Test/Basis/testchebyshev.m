classdef testchebyshev < matlab.unittest.TestCase
    properties (TestParameter)
        n = {1,2,3};
        nMonomials = {1,2,3};
        nShifts = {1,2,3};
    end

    methods (Test)
        % Test the method:
        function shifttest(testCase)
            A = [1 2 3];
            B = [1 1 1; 2 0 2];
            actual = chebyshev(A,B);
            expected1 = [2 3 4; 3 2 5; 0 3 4; 1 2 5; 2 1 4; 3 2 5; ...
                0 1 4; 1 2 5; 2 3 2; 3 2 1; 0 3 2; 1 2 1; 2 1 2; ...
                3 2 1; 0 1 2; 1 2 1];
            expected2 = (1/8)*ones(16,1);
            expected = [expected2 expected1];
            testCase.verifyEqual(actual,expected);
        end
            
        function recursiontest(testCase)
            [~,alpha,beta,phi0,phi1] = chebyshev(ones(1,3),ones(1,3));
            expectedAlpha = [2 1 1 1];
            expectedBeta = [-1 0 0 0];
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
            actual = chebyshev(A,B);
            dilation = 2^n;
            expectedSize = [dilation*nShifts*nMonomials, n+1];
            testCase.verifySize(actual,expectedSize);
        end
    end
end
