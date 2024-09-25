classdef testcontraction < matlab.unittest.TestCase
    methods (Test)
        % Test the method:
        function contractiontest(testCase)
            p = contraction([6 0 4 3 2 1],2,2);
            testCase.verifyEqual(p,[1 0 0; 2 0 1; 3 1 0; 4 0 2; 6 2 0]);
        end

        function tolerancetest(testCase)
            p = contraction([6 0 1e-6 3 2 1],2,2,1e-5);
            testCase.verifyEqual(p,[1 0 0; 2 0 1; 3 1 0; 6 2 0]);
        end

        function complextest(testCase)
            p = contraction([1 2i 3-2i 1e-6 0 0],2,2,1e-5);
            actual = sortrows(p);
            expected = sortrows([1 2 0; 2i 1 1; 3-2i 0 2]);
            testCase.verifyEqual(actual,expected);
        end

        function ordertest(testCase)
            p = contraction([1 2 3 4 0 6],2,2,@grnlex);
            actual = sortrows(p);
            expected = sortrows([1 0 2; 2 1 1; 3 2 0; 4 0 1; 6 0 0]);
            testCase.verifyEqual(actual,expected);
        end

        function combinationtest(testCase)
            K = monomials(3,3);
            coef = 1:length(K);
            A = [coef' K];
            B = circshift(A,3);
            p = contraction(expansion(B),3,3);
            testCase.verifyEqual(p,A);
        end
    end
end
