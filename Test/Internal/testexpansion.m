classdef testexpansion < matlab.unittest.TestCase
    methods (Test)
        % Test the method:
        function univariatetest(testCase)
            K = monomials(3,3);
            coef = 1:length(K);
            A = [coef' K];
            B = circshift(A,3);
            testCase.verifyTrue(all(expansion(B) == fliplr(coef)));
        end

        function tolerancetest(testCase)
            P = [1 0 0; 1e-5 0 1; 1 1 0; 2 0 3];
            actual = expansion(P,1e-4);
            expected = [0 0 0 2 0 0 0 1 0 1];
            testCase.verifyEqual(actual,expected);
        end

        function ordertest(testCase)
            P = [1 2 0; 2 0 0; 3 0 1; 4 1 0];
            actual = expansion(P,@grnlex);
            expected = [0 0 1 3 4 2];
            testCase.verifyEqual(actual,expected);
        end

        function complextest(testCase)
            P = [1i 0 0; 1+2i 0 1; 1-3i 1 0; 2 1 1];
            actual = expansion(P);
            expected = [0 2 0 1-3i 1+2i 1i];
            testCase.verifyEqual(actual,expected);
        end
    end
end
