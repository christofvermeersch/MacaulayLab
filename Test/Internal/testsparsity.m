classdef testsparsity < matlab.unittest.TestCase
    methods (Test)
        % Test the method:
        function sparsetest(testCase)
            p = [1 0 2; 1 1 1];
            q = [1 2 0];
            system = systemstruct({p,q});
            actual = sparsity(system);
            expected = 0.25;
            testCase.verifyEqual(actual,expected);
        end

        function fulltest(testCase)
            mep = randmep(2,2,5,4);
            actual = sparsity(mep);
            expected = 1.0;
            testCase.verifyEqual(actual,expected);
        end

        function tolerancetest(testCase)
            A00 = [1 2; 3 4];
            A10 = [1 0.01i; -0.01i 0];
            A01 = [1e-5 1; 2 3];
            mep = mepstruct({A00,A10,A01},[0 0; 1 0; 0 1]);
            actual1 = sparsity(mep);
            actual2 = sparsity(mep,1e-16);
            actual3 = sparsity(mep,1e-4);
            expected1 = 11/12;
            expected2 = 11/12;
            expected3 = 10/12;
            testCase.verifyEqual(actual1,expected1);
            testCase.verifyEqual(actual2,expected2);
            testCase.verifyEqual(actual3,expected3);
        end
    end
end
