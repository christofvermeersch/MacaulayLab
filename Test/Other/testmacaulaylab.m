classdef testmacaulaylab < matlab.unittest.TestCase
    properties
        tolerance = 1e-10
    end

    methods (Test)
        % Test the method:
        function systemtest(testCase)
            p = [-1 2 0; 2 1 1; 1 0 2; 5 1 0; -3 0 1; -4 0 0];
            q = [1 2 0; 2 1 1; 1 0 2; -1 0 0];
            system = systemstruct({p,q});
            actual = sortrows(num(macaulaylab(system)));
            expected = [0 -1; 1 0; 3 -2; 4 -5];
            testCase.verifyLessThanOrEqual(actual-expected,testCase.tolerance);
        end

        function peptest(testCase)
            A2 = diag([3 1 3 1]);
            A1 = [0.4 0 -0.3 0; 0 0 0 0; -0.3 0 0.5 -0.2; 0 0 -0.2 0.2];
            A0 = [-7 2 4  0; 2 -4 2 0; 4 2 -9 3; 0 0 3 -3];
            system = mepstruct({A0,A1,A2},[0; 1 ; 2]);
            actual = sort(num(macaulaylab(system)));
            expected = sort(polyeig(A0,A1,A2));
            testCase.verifyLessThanOrEqual(actual-expected,testCase.tolerance);
        end
    end
end
