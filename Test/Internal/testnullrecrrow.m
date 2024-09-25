classdef testnullrecrrow < matlab.unittest.TestCase
    properties (Constant)
        tolerance = 1e-10;
    end

    methods (Test)
        % Test the method:
        function standardtest(testCase)
            A = randn(10,5)*randn(5,10);
            B = randn(10,3)*randn(3,10);
            C = randn(10,2)*randn(2,10);
            
            Z = nullrecrrow(null([A; B]), C);
            testCase.verifySize([A; B; C]*Z,[30 0]);
        end
        
        function exceptionaltest(testCase)
            A = randn(10,5)*randn(5,10);
            B = randn(10,3)*randn(3,10);
            C = randn(10,2)*randn(2,10);
            
            Z = nullrecrrow(null([A; B; C]), A);
            testCase.verifySize([A; B; C; A]*Z, [40 0]);
        end
    end
end