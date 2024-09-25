classdef testresiduals < matlab.unittest.TestCase
    properties (Constant)
        tolerance = 1e-10;
    end

    methods (Test)
        % Test the method:
        function polynomialtest(testCase)
            P = cell(2,1);
            P{1} = [-1 2 0; 2 1 1; 1 0 2; 5 1 0; -3 0 1; -4 0 0];
            P{2} = [1 2 0; 2 1 1; 1 0 2; -1 0 0];
            system = systemstruct(P);
            
            [maxres, res] = residuals([1 0; 1 1; -1 2; 0 0],system);
            testCase.verifyEqual(maxres,16);
            testCase.verifyEqual(res,[0; 3; 16; norm([-4 -1])]);
        end

        function meptest(testCase)
            P = cell(3,1);
            P{1} = [1 1; 1 1; 1 1];
            P{2} = [1 0; 1 0; 1 0];
            P{3} = [0 2; 0 2; 0 2];
            x = [2 1];
            maxres = residuals(x,mepstruct(P,[0 0; 0 1; 1 0]));
            testCase.verifyLessThanOrEqual(maxres,testCase.tolerance);
        end
    end
end
