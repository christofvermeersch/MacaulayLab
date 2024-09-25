classdef teststdmonomials < matlab.unittest.TestCase
    methods (Test)
        % Test the method:
        function exacttest(testCase)
            A = randn(3,2)*randn(2,5);
            B = randn(9,2)*randn(2,5);
            C = randn(30,3)*randn(3,5);
            Z = [A; B; A + 10e-6*randn(3,5); B; A ; A; C];
            c = stdmonomials(Z);
            testCase.verifyEqual(c,[1; 2; 4; 5; 13]);
        end

        function tolerancetest(testCase)
            A = randn(3,2)*randn(2,5);
            B = randn(9,2)*randn(2,5);
            C = randn(30,3)*randn(3,5);
            Z = [A; B; A + 10e-6*randn(3,5); B; A ; A; C];
            c = stdmonomials(Z,10e-5);
            testCase.verifyEqual(c,[1; 2; 4; 5; 31]);
        end
    end
end