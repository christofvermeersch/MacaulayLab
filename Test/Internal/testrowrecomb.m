classdef testrowrecomb < matlab.unittest.TestCase
    methods (Test)
        % Test the method:
        function simpetest(testCase)
            actual = rowrecomb(2,2,[1 2],[2 0 1; 3 1 0]);
            expected = [1 0 0 0 0 0; 0 1 0 0 0 0; 0 2 3 0 0 0; ...
                0 0 0 2 3 0; 0 0 0 0 1 0; 0 0 0 0 0 1];
            testCase.verifyEqual(full(actual),expected);
        end

        function blocksizetest(testCase)
            actual = rowrecomb(2,2,[1 2],[1 0 1],2);
            expected = eye(12);
            testCase.verifyEqual(full(actual),expected);
        end
        
        function ordertest(testCase)
            actual = rowrecomb(2,2,[1 2],[2 1 0; 3 0 1],@grnlex);
            expected = [1 0 0 0 0 0; 0 1 0 0 0 0; 0 2 3 0 0 0; ...
                0 0 0 2 3 0; 0 0 0 0 1 0; 0 0 0 0 0 1];
            testCase.verifyEqual(full(actual),expected);
        end

        function basistest(testCase)
            actual = rowrecomb(4,2,[1,5],[1 1 1],@chebyshev);
            expected = zeros(15,15);
            expected(1,1) = 1;
            expected(2,5) = 1;
            expected(3,1) = 0.25;
            expected(3,4) = 0.25;
            expected(3,6) = 0.25;
            expected(3,13) = 0.25;
            expected(4,2) = 1;
            expected(5,3) = 1;
            expected(6,6) = 1;
            expected(7,7) = 1;
            expected(8,8) = 1;
            expected(9,9) = 1;
            expected(10,10) = 1;
            expected(11,11) = 1;
            expected(12,12) = 1;
            expected(13,13) = 1;
            expected(14,14) = 1;
            expected(15,15) = 1;
            testCase.verifyEqual(full(actual),expected);
        end

        function complextest(testCase)
            actual = rowrecomb(2,2,[1 2],[2+3i 0 1; -2i 1 0]);
            expected = [1 0 0 0 0 0; 0 1 0 0 0 0; 0 2+3i -2i 0 0 0; ...
                0 0 0 2+3i -2i 0; 0 0 0 0 1 0; 0 0 0 0 0 1];
            testCase.verifyEqual(full(actual),expected);
        end
    end
end
