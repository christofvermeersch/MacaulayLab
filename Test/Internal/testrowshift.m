classdef testrowshift < matlab.unittest.TestCase
    methods (Test)
        % Test the method:
        function simpletest(testCase)
            actual = rowshift(2,2,[1 2],[3 1 0; 2 0 1]);
            expected = [0 2 3 0 0 0; 0 0 0 2 3 0];
            testCase.verifyEqual(full(actual),expected);
        end

        function blocksizetest(testCase)
            actual = rowshift(2,2,[1 2],[1 0 1],2);
            expected = [0 0 1 0 0 0 0 0 0 0 0 0; ...
                0 0 0 1 0 0 0 0 0 0 0 0];
            testCase.verifyEqual(full(actual),expected);
        end
        
        function ordertest(testCase)
            actual = rowshift(2,2,[1 2],[2 1 0; 3 0 1],@grnlex);
            expected = [0 2 3 0 0 0; 0 0 0 2 3 0];
            testCase.verifyEqual(full(actual),expected);
        end

        function basistest(testCase)
            actual = rowshift(4,2,5,[1 1 1],@chebyshev);
            expected = zeros(1,15);
            expected(1,1) = 0.25;
            expected(1,13) = 0.25;
            expected(1,4) = 0.25;
            expected(1,6) = 0.25;
            testCase.verifyEqual(full(actual),expected);
        end

        function complextest(testCase)
            actual = rowshift(2,2,[1 2],[2+3i 1 0; -2i 0 1]);
            expected = [0 -2i 2+3i 0 0 0; 0 0 0 -2i 2+3i 0];
            testCase.verifyEqual(full(actual),expected);
        end

        function thesistest(testCase)
            V = vandermonde(2,2,[2 -1; 4 1]);
            rows = [1 2];
            S1 = rowshift(2,2,rows,[1 0 0]);
            Sx1 = rowshift(2,2,rows,[1 1 0]);
            Sg = rowshift(2,2,rows,[2 1 0; 3 0 1]);
            left = S1*V*diag([2 4]);
            right = Sx1*V;
            testCase.verifyEqual(left,right);
            left = S1*V*diag([1 11]);
            right = Sg*V;
            testCase.verifyEqual(left,right);
        end
    end
end
