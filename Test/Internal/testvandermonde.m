classdef testvandermonde < matlab.unittest.TestCase
    methods (Test)
        % Test the method:
        function thesistest(testCase)
            actual = vandermonde(2,2,[2 -1; 4 1]);
            expected = [1 1; -1 1; 2 4; 1 1; -2 4; 4 16];
            testCase.verifyEqual(actual,expected);
        end

        function nonblockedtest(testCase)
            actual = vandermonde(3,2,[-2 3; 4 5]);
            expected = [1 1; 3 5; -2 4; 9 25; -6 20; 4 16; 27 125; ...
                -18 100; 12 80; -8 64];
            testCase.verifyEqual(actual,expected);
        end

        function blockedtest(testCase)
            actual = vandermonde(2,2,[2 1],[1 2]);
            expected = [2; 1; 2; 1; 4; 2; 2; 1; 4; 2; 8; 4]; 
            testCase.verifyEqual(actual,expected);
        end

        function ordertest(testCase)
            actual1 = vandermonde(2,2,[1 2],[1 2],@grnlex);
            expected1 = [1; 2; 1; 2; 2; 4; 1; 2; 2; 4; 4; 8]; 
            testCase.verifyEqual(actual1,expected1);
            actual2 = vandermonde(3,2,[-2 3; 4 5],@grnlex);
            expected2 = [1 1; -2 4; 3 5; 4 16; -6 20; 9 25; -8 64; ...
                12 80; -18 100; 27 125];
            testCase.verifyEqual(actual2,expected2);
        end

        function complextest(testCase)
            actual = vandermonde(2,3,[1i 1+1i 2-2i]);
            expected = [1; 2-2i; 1+1i; 1i; -8i; 4; 2+2i; 2i; -1+1i; -1];
            testCase.verifyEqual(actual,expected);
        end
    end
end
