classdef testmonomials < matlab.unittest.TestCase
    methods (Test)
        % Test the method:
        function matrix1test(testCase)
            testCase.verifyEqual(monomials(0,4),[0 0 0 0]);
        end

        function matrix2test(testCase)
            testCase.verifyEqual(monomials(2,2), ...
                [0 0; 0 1; 1 0; 0 2; 1 1; 2 0]);
        end

        function matrix3test(testCase)
            testCase.verifyEqual(monomials(1,3), ...
                [0 0 0; 0 0 1; 0 1 0; 1 0 0]);
        end

        function ordertest(testCase)
            M1 = monomials(2,3,@grevlex);
            M2 = monomials(2,3,@grlex);
            M3 = monomials(2,3,@grnlex);
            M4 = monomials(2,3,@grvlex);
            N1 = [0 0 0; 0 0 1; 0 1 0; 1 0 0; 0 0 2; 0 1 1; 1 0 1; ...
                0 2 0; 1 1 0; 2 0 0];
            N2 = [0 0 0; 0 0 1; 0 1 0; 1 0 0; 0 0 2; 0 1 1; 0 2 0; ...
                1 0 1; 1 1 0; 2 0 0];
            N3 = [0 0 0; 1 0 0; 0 1 0; 0 0 1; 2 0 0; 1 1 0; 1 0 1; ...
               0 2 0; 0 1 1; 0 0 2];
            N4 = [0 0 0; 1 0 0; 0 1 0; 0 0 1; 2 0 0; 1 1 0; 0 2 0; ...
                1 0 1; 0 1 1; 0 0 2];
            testCase.verifyEqual(M1,N1); 
            testCase.verifyEqual(M2,N2); 
            testCase.verifyEqual(M3,N3); 
            testCase.verifyEqual(M4,N4); 
        end

        function blockedtest(testCase)
            M = monomials(3,2,2);
            N = [0 0; 0 0; 0 1; 0 1; 1 0; 1 0; 0 2; 0 2; 1 1; 1 1; 2 0; ...
                2 0; 0 3; 0 3; 1 2; 1 2; 2 1; 2 1; 3 0; 3 0];
            testCase.verifyEqual(M,N);
        end
    end
end
