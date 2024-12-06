classdef testsystemstruct < matlab.unittest.TestCase
    methods (Test)
        % Test the class constructor:
        function constructortest(testCase)
            p = [1 5 0; 2 2 2];
            q = [3 4 2; 4 1 1];
            system = systemstruct({p,q});
            testCase.verifyClass(system,"systemstruct");
            testCase.verifyInstanceOf(system,"problemstruct");
            testCase.verifyEqual(system.eqs,{p,q});
            testCase.verifyEqual(system.coef,{[1; 2];[3; 4]});
            testCase.verifyEqual(system.supp,{[5 0; 2 2];[4 2; 1 1]});
            testCase.verifyEqual(system.basis,"unknown");
        end

        function symbolictest(testCase)
            error("Test is not yet implemented.")
        end

        function maptest(testCase)
            error("Test is not yet implemented.")
        end
    end
end
